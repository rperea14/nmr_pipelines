classdef dwi_HAB < dwiMRI_Session
    %%  classdef dwi_ADRC < dwiMRI_Session
    %%  This class is a subclass of its parent class dwi_MRI_Session.m
    %%  (where it will inherent other methods).
    %%  Created by: Rodrigo Perea rpereacamargo@mgh.harvard.edu
    %%              Aaron Schultz aschultz@martinos.org
    %%              
    %%
    %%      Dependencies (and tested in):
    %%          -FreeSurfer v6.0
    %%          -SPM8
    %%          -Ants tools 1.7.1
    %%          -DSI_studio_vApr19_2017
    %%          -FSL 5.0.9
    %%  *Only filesep wit '/' are used in the properties class declaration.
    %%   Besides these, all should be any operating system compatible (tested in CentOS Linux)
    properties
        %Properties initialized in parent class dwiMRI_Session() and 
        %And assgiend values in method: obj.setMyParams() which is called
        %in the class constructor dwi_HAB().
    end
    
    methods
        function obj = dwi_HAB(sessionname,opt)
            %%%%%%%%%%%%%%%%%%%%%%%%
            %Initialize Parameters:
            obj.setMyParams();
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %Initialize root variables:
            obj.sessionname = sessionname;
            obj.root = [obj.session_dir sessionname '/DWIs/'];
            obj.dcm_location= [ obj.dcm_location sessionname filesep ];
            obj.session_location= [ obj.session_dir sessionname filesep ];
            %Assign Project ID:
            obj.projectID='HAB';
            
            %Check if you are in the right project (if not, probably
            %obj.session_location doesn't exist:
            if ~exist([obj.session_location ],'dir')
                error(['Sesion directory: ' obj.session_location ' doesn''t exist. Are you sure you are in the right Project ID?']);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %If the folder <XX>/DWIs/ does not exist, then create it!
            if exist(obj.root,'dir')==0
                obj.make_root();
            end
            obj.objectHome = obj.root ;
            if exist([obj.objectHome filesep sessionname '.mat'],'file') > 0
                load([obj.objectHome filesep sessionname '.mat']);
                oldroot = obj.root;
                obj.wasLoaded = true;
            else
                donothing='';
                %obj.setMyParams;
            end
            
            
            %Reinitialize params if obj.FSL_dir is empty:
            if isempty(obj.FSL_dir)
               obj.setMyParams(); 
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %CHECK CHANGES MADE FROM DWIs_XX/Sessions/DWIs to
            %Sessions/DWIs:
            if strcmp(strtrim(['/cluster/sperling/HAB/Project1/DWIs_30b700/Sessions/' obj.sessionname '/DWIs/' ] ),obj.root)
                display('Changing sessions folder...');
                %Then apply the mod to cluster...
                newobj = replaceObjText_v2(obj,{obj.root},{['/cluster/sperling/HAB/Project1/Sessions/' obj.sessionname '/DWIs/']});
                obj=newobj;
                obj.objectHome = obj.root ;
                display(['CHANGING: ''<HAB1>/DWIs_30b700/Sessions/' obj.sessionname  'DWIs'' was changed to ''<HAB1>/Sessions/' obj.sessionname   '/DWIs/'' ']);
                obj.resave();
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %CHECK CHANGES MADE FROM /eris to /cluster
            if strcmp(strtrim(obj.FS_location),'/eris/bang/HAB_Project1/FreeSurferv6.0')
                display('Changing eris to cluster folder...');
                %Then apply the mod to cluster...
                newobj = replaceObjText_v2(obj,{'eris'},{'cluster'});
                obj=newobj;
                obj.objectHome = obj.root ;
                display('CHANGING: ''eris'' was changed to ''cluster'' ');
                obj.resave();
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Check if *.nii.gz files exist, if not get them from DCM2nii:
            obj.rawfiles = dir_wfp([obj.root 'Orig/*.nii.gz' ] );
            if isempty(obj.rawfiles)
                RunFlag=true;
                obj.getDCM2nii(RunFlag);
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Continue with CommonPreProc
            obj.CommonPreProc();
            %Continue with CommonPostProc
            obj.CommonPostProc();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isdeployed()
                display('Diffusion Pre- and Post- processing completed')
                fprintf('\n\n *If running the compiled version and for completion, \n make sure you run (from MATLAB terminal):')
                fprintf('\n * 1. Load the object ~~> load(<Sessions>/<SUBJID>/DWIs/<SUBJID>.mat) \n')
                fprintf('\n * 2. Run             ~~> obj.proc_T1toDWI() \n')
                fprintf('\n * 3. Run             ~~> obj.trkland_fx() \n')
            end
            
        end
       
        function obj = CommonPreProc(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Get DCM2NII File location:
            RunFlag=false; %here, we only allocate variable, not run proc_dcm2nii
            obj.getDCM2nii(RunFlag);
            %Selecting the raw files
            obj.rawfiles = strtrim(dir_wfp([obj.root 'Orig/*.nii.gz' ] ));
      
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For BET2:
            obj.Params.Bet2.in.movefiles = ['..' filesep '01_Bet'];
            obj.Params.Bet2.in.fracthrsh = 0.4;
            obj.Params.Bet2.in.fn = obj.rawfiles;
            
            obj.proc_bet2();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For EDDY:
            obj.Params.Eddy.in.movefiles = ['..' filesep '02_Eddy'];
            obj.Params.Eddy.in.fn=obj.rawfiles;
            obj.Params.Eddy.in.bvals=strrep(obj.rawfiles,'.nii.gz','.bvals');
            obj.Params.Eddy.in.bvecs=strrep(obj.rawfiles,'.nii.gz','.voxel_space.bvecs');
            obj.Params.Eddy.in.mask = obj.Params.Bet2.out.mask;
            obj.Params.Eddy.in.index= ones(1,35) ; %for 35 volumes
            obj.Params.Eddy.in.acqp= [ 0 -1 0 0.102]; %based on HAB diff sequence
            
            obj.proc_eddy();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For B0mean:
            obj.Params.B0mean.in.movefiles = ['..' filesep '03_B0mean'];
            obj.Params.B0mean.in.fn=obj.Params.Eddy.out.fn;
            obj.Params.B0mean.in.b0_nvols=5;
            
            obj.proc_meanb0();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %To calculate mean motion based on edddy output:
            obj.Params.EddyMotion.in.movefiles = ['..' filesep '03_EddyMotion']; 
            obj.Params.EddyMotion.in.fn_eddy = obj.Params.Eddy.out.fn ;
            
            obj.proc_get_eddymotion();
                        
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %To generate a mask after Eddy:
            %This step will 1) define a better mask if eddy affected the
            %movement of the head and 2) remove issues known to happen at
            %the edges of the brain when using the --wls option in dtifit!\
            %Reference: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;6eb4d787.1610
            obj.Params.MaskAfterEddy.in.movefiles = ['..' filesep '04_MaskAfterEddy'];
            obj.Params.MaskAfterEddy.in.fn=obj.Params.B0mean.out.fn;
            obj.Params.MaskAfterEddy.in.prefix = 'after_eddy';
            obj.Params.MaskAfterEddy.in.fracthrsh = '0.4';
            
            obj.proc_mask_after_eddy();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For DTIFIT:
            %We will use the --wls option as it seems to improve the fit of
            %the diffusion tensor model and negativity values due to noise
            %REF: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;735f6320.1409
            obj.Params.Dtifit.in.movefiles = [ '..' filesep '05_Recon_dtifit' ];
            obj.Params.Dtifit.in.fn = obj.Params. Eddy.out.fn;
            obj.Params.Dtifit.in.prefix = 'DTIFIT_FSLv509' ; %Double check this so you prefix the version of FSL!
            obj.Params.Dtifit.in.bvecs = obj.Params.Eddy.out.bvecs;
            obj.Params.Dtifit.in.bvals = obj.Params.Eddy.in.bvals;
            obj.Params.Dtifit.in.mask = obj.Params.MaskAfterEddy.out.finalmask;
            
            obj.proc_dtifit();
         
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For GQI:
            obj.Params.GQI.in.movefiles = [ '..' filesep '05_Recon_gqi' ];
            obj.Params.GQI.in.fn = obj.Params. Eddy.out.fn;
            obj.Params.GQI.in.mask = obj.Params.MaskAfterEddy.out.finalmask;
            obj.Params.GQI.in.bvecs = obj.Params.Eddy.out.bvecs;
            obj.Params.GQI.in.bvals = obj.Params.Eddy.in.bvals;
            
            obj.Params.GQI.in.prefix = obj.dsistudio_version ; %Double check this so you prefix the version of DSISTUDIO!
            obj.Params.GQI.out.export = 'gfa,nqa0,nqa1';
            
            obj.Params.GQI.in.method = '4';    %for gqi model
            obj.Params.GQI.in.num_fiber = '3'; %modeling 3 fiber population
            obj.Params.GQI.in.param0 = '1.25'; %default parameter for gqi
            
            obj.proc_gqi();
      
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For AntsReg:
            obj.Params.AntsReg.in.movefiles = ['..' filesep '06_Ants_CoReg' ];
            obj.Params.AntsReg.in.run_ants = '/cluster/bang/sw/ANTs/ANTS-2.1.0-redhat';
            obj.Params.AntsReg.in.prefix = 'Antsv201_2_HABn272_';
            obj.Params.AntsReg.in.fn = obj.Params.Dtifit.out.FA;
            obj.Params.AntsReg.in.ref = obj.HABn272_meanFA;
            obj.Params.AntsReg.in.antsparams = ' -d 3 -n 4 -t s -r 4 -p d ' ;
            
            obj.proc_antsreg(); 
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For Skeletonize:
            obj.Params.Skeletonize.in.movefiles = ['..' filesep '07_Skeletonize' ];
            obj.Params.Skeletonize.in.fn = obj.Params.AntsReg.out.fn;
            obj.Params.Skeletonize.in.meanFA = obj.HABn272_meanFA;
            obj.Params.Skeletonize.in.skel_dst = obj.HABn272_meanFA_skel_dst;
            obj.Params.Skeletonize.in.thr = '0.3';
            obj.Params.Skeletonize.in.ref_region = obj.ref_region;
            obj.Params.Skeletonize.in.prefix = [ 'FSLv509_skelHABn272' ]; %check fsl versio you are calling!
            
            obj.proc_skeletonize();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For Skel_TOI:
            obj.Params.Skel_TOI.in.location = obj.skeltoi_location;
            obj.Params.Skel_TOI.in.masks = obj.skeltoi_tois ;
            obj.Params.Skel_TOI.in.ext = '.nii.gz' ;
            obj.Params.Skel_TOI.in.suffix = '_n272TMP';
           
            obj.proc_getskeltois();
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Upload skeletons and diffusion motion to DataCentral:
            obj.UploadData_DWI(); %*Make sure you have a connDB.p in your path for this to work - rdp 4/12/18!
    
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For FreeSurfer recon-all
            [ ~, tmp_shell ] = system('echo $SHELL ');
            [~, tmpb ] = fileparts(tmp_shell);
            obj.Params.FreeSurfer.shell = strtrim(tmpb); %strtrim(tmpshell);
            
            obj.Params.FreeSurfer.dir = obj.FS_location ;
            obj.Params.FreeSurfer.init_location = obj.init_FS; 
            %Retrieving a T1 scan:
            [sys_error, obj.Params.FreeSurfer.in.T1raw ] = system(['ls ' obj.session_location 'MPRAGE' filesep '*.mgz | head -1' ]);
            if sys_error ~= 0 %No problem, we get the T1 the continue...
                fprintf(['\nError when finding the T1:'  obj.Params.FreeSurfer.in.T1raw  '\n'])
            end
            %Retrieving a T2 scan:
            [sys_error, obj.Params.FreeSurfer.in.T2raw ] = system(['ls ' obj.session_location 'other' filesep '*T2* | head -1' ]);
            if sys_error ~= 0 %No problem, we get the T1 the continue...
                fprintf(['\nNo T2 found:'  obj.Params.FreeSurfer.in.T2raw  '\n'])
                obj.Params.FreeSurfer.in.T2exist=false;
            else
                obj.Params.FreeSurfer.in.T2exist=true;
            end
            obj.Params.FreeSurfer.out.aparcaseg = [ obj.Params.FreeSurfer.dir ...
                filesep obj.sessionname filesep 'mri' filesep 'aparc+aseg.mgz' ] ;
            
            obj.proc_getFreeSurfer();
             
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %FS2dwi (move aparc and aparc2009 segmes to dwi space):
            obj.Params.FS2dwi.in.movefiles = ['..' filesep '05_FS2dwi' ];
            obj.Params.FS2dwi.in.b0 = obj.Params.B0mean.out.fn ; 
            obj.Params.FS2dwi.in.aparcaseg = obj.Params.FreeSurfer.out.aparcaseg ; 
            obj.Params.FS2dwi.in.tmpfile_aparcaseg = [ obj.dependencies_dir  filesep 'FS_DEPS' filesep  'FS_aparc.txt' ] ; 
            obj.Params.FS2dwi.in.tmpfile_aparcaseg2009 = [ obj.dependencies_dir  filesep 'FS_DEPS' filesep  'FS_aparc2009.txt' ] ; 
            obj.Params.FS2dwi.in.tmpfile_hippo_bil = [ obj.dependencies_dir  filesep 'FS_DEPS' filesep  'FS_hippolabels_bil.txt' ] ;
            obj.Params.FS2dwi.in.aparcaseg2009 = ...
                strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','aparc.a2009s+aseg')); 
            %*A possible error is the naming convention when only a T1 was used!
            %Using this will allow us to automatically select hippo-fields without discriminating whetehr they come from only a T1 or using T1-T2:            
            [~ , tmp_hippo_L ] = system(['ls ' strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','lh.hippoSfLabels*') ' | tail -1 ']);
            obj.Params.FS2dwi.in.hippofield_left = strtrim(tmp_hippo_L);
            [~ , tmp_hippo_R ] = system(['ls ' strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','rh.hippoSfLabels*') ' | tail -1 ']);
            obj.Params.FS2dwi.in.hippofield_right = strtrim(tmp_hippo_R);
            clear tmp_hippo_L tmp_hippo_R;
%             obj.Params.FS2dwi.in.hippofield_left = ...
%                 strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','lh.hippoSfLabels-T1-T2.v10.FSvoxelSpace')); 
%             obj.Params.FS2dwi.in.hippofield_right = ...
%                 strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','rh.hippoSfLabels-T1-T2.v10.FSvoxelSpace')); 
%             
            obj.proc_FS2dwi();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Extracting the value sfrom FreeSurfer:
            obj.getdata_FreeSurfer();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Coregistering the T1 to B0:
            obj.Params.T1toDWI.in.movefiles = './05_T1toDWI/';
            obj.Params.T1toDWI.in.b0 = strtrim(obj.Params.B0mean.out.fn{1}); %strtrim will remove the leading \n special character
            obj.Params.T1toDWI.in.T1 = strtrim(obj.Params.FreeSurfer.in.T1);
            
            %obj.proc_T1toDWI(); %!!*NOT IN COMPILED VERSION DUE TO ERROR
            %READING PRIVATE FUNCTIONS (read_hdr(), part of spm12)
            
            %~~~Continue with obj.CommonPostProc...>
        end
        
        function obj = CommonPostProc(obj) 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %TRKLAND
            obj.Trkland.root = [ obj.root  'post_TRKLAND' filesep ];
            %FX_TRKLAND:
            for tohide=1:1
                %b0 params:
                obj.Trkland.fx.in.b0 = obj.Params.B0mean.out.fn{end}; % obj.Params.MaskAfterEddy.out.brainonly; % obj.Params.B0mean.out.fn{end}; --> BAD 
                obj.Trkland.fx.in.movefiles = ['..' filesep 'post_TRKLAND' ];
               
                %Based on orientation, we will chooose a specific template
                %(usually RAS) 
                [~, Ori ] = system(['mri_info ' obj.Trkland.fx.in.b0 ' | grep Orientation | awk ''{print $3}'''] );
                obj.Trkland.fx.tmp.ori = strtrim(Ori);
                clear Ori
                if strcmp(obj.Trkland.fx.tmp.ori,'LPS') || strcmp(obj.Trkland.fx.tmp.ori,'LAS')
                    obj.Trkland.fx.tmp.b0 = [ obj.fx_template_dir 'LPS_141210_8CS00178_b0.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roa_solid_bil =[ obj.fx_template_dir 'LPS_TMP_178_bil_fx_dil11.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roa_solid_lh = [ obj.fx_template_dir 'LPS_TMP_178_lh_fx_dil11.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roa_solid_rh = [ obj.fx_template_dir 'LPS_TMP_178_rh_fx_dil11.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roi_bil = [ obj.fx_template_dir 'LPS_TMP_178_bil_fx_dil.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roi_lh = [ obj.fx_template_dir 'LPS_TMP_178_lh_fx_dil.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roi_rh = [ obj.fx_template_dir 'LPS_TMP_178_rh_fx_dil.nii.gz' ] ;
                elseif strcmp(obj.Trkland.fx.tmp.ori,'RAS')
                    obj.Trkland.fx.tmp.b0 = [ obj.fx_template_dir 'RAS_141210_8CS00178_b0.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roa_solid_bil =[ obj.fx_template_dir 'RAS_TMP_178_bil_fx_dil11.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roa_solid_lh = [ obj.fx_template_dir 'RAS_TMP_178_lh_fx_dil11.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roa_solid_rh = [ obj.fx_template_dir 'RAS_TMP_178_rh_fx_dil11.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roi_bil = [ obj.fx_template_dir 'RAS_TMP_178_bil_fx_dil.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roi_lh = [ obj.fx_template_dir 'RAS_TMP_178_lh_fx_dil.nii.gz' ] ;
                    obj.Trkland.fx.tmp.roi_rh = [ obj.fx_template_dir 'RAS_TMP_178_rh_fx_dil.nii.gz' ] ;
                else
                    error('In trkland_fx() Init: Cannot find the right fornix template orientation to use for co-registration. Quitting...');
                end
                %IN PARAMS:
                %Hippocampi:
                obj.Trkland.fx.in.hippo_lh =  strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Left-Hippocampus.nii.gz');
                obj.Trkland.fx.in.hippo_rh =  strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Right-Hippocampus.nii.gz');
                %Thalami:
                obj.Trkland.fx.in.thalamus_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Left-Thalamus-Proper.nii.gz');
                obj.Trkland.fx.in.thalamus_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Right-Thalamus-Proper.nii.gz');
                %tmp2b0s params:
                obj.Trkland.fx.in.fn_tmp2b0 =  [ obj.Trkland.root 'fx_tmp2b0.nii.gz' ];
                obj.Trkland.fx.in.tmp2b0_matfile = [ obj.Trkland.root 'fx_tmp2b0.mat'];
                %bil params:
                obj.Trkland.fx.in.roi_bil = [ obj.Trkland.root 'fx_roi_bil.nii.gz'];
                obj.Trkland.fx.in.roa_bil_solid = [ obj.Trkland.root 'fx_roa_bil_solid.nii.gz'];
                obj.Trkland.fx.in.roa_bil_ero =  [ obj.Trkland.root 'fx_roa_bil_ero.nii.gz'];
                obj.Trkland.fx.in.roa_bil_hollow = [ obj.Trkland.root 'fx_roa_bil_hollow.nii.gz'];
                %lh params:
                obj.Trkland.fx.in.roi_lh_hippo = strrep(obj.Params.FS2dwi.out.fn_aparc2009, ...
                    'dwi_aparc.a2009+aseg',[ 'aparc2009_aseg' filesep 'dwi_fs_Left-Hippocampus' ]);
                obj.Trkland.fx.in.roi_lh = [ obj.Trkland.root 'fx_roi_lh.nii.gz'];
                obj.Trkland.fx.in.roa_lh_solid = [ obj.Trkland.root 'fx_roa_lh_solid.nii.gz'];
                obj.Trkland.fx.in.roa_lh_ero = [ obj.Trkland.root 'fx_roa_lh_ero.nii.gz'];
                obj.Trkland.fx.in.roa_lh_hollow = [ obj.Trkland.root 'fx_roa_lh_hollow.nii.gz'];
                %rh params:
                obj.Trkland.fx.in.roi_rh_hippo = strrep(obj.Params.FS2dwi.out.fn_aparc2009, ...
                    'dwi_aparc.a2009+aseg',[ 'aparc2009_aseg' filesep 'dwi_fs_Right-Hippocampus' ]);
                obj.Trkland.fx.in.roi_rh = [ obj.Trkland.root 'fx_roi_rh.nii.gz'];
                obj.Trkland.fx.in.roa_rh_solid = [ obj.Trkland.root 'fx_roa_rh_solid.nii.gz'];
                obj.Trkland.fx.in.roa_rh_ero = [ obj.Trkland.root 'fx_roa_rh_ero.nii.gz'];
                obj.Trkland.fx.in.roa_rh_hollow = [ obj.Trkland.root 'fx_roa_rh_hollow.nii.gz'];
                %fib params:
                obj.Trkland.fx.in.fib =strtrim(obj.Params.GQI.out.fibs_fn{end});
                if exist(obj.Trkland.fx.in.fib) == 0 ; error('No fib found in variable: trkland.trks.fx.in.fib. Please check!') ; end
                
                %Interpolation n:
                obj.Trkland.fx.in.n_interp=40; %According to average value on previous studies in connectome!
                %obj.trkland_fx(); %!!*NOT IN COMPILED VERSION DUE TO ERROR
            %READING PRIVATE FUNCTIONS (read_hdr(), part of spm12)
               
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %TRACULA (and implicit functionality of bedpostx):
            for tohide=1:1
                obj.Params.Tracula.in.movefiles = ['..' filesep 'post_TRACULA' ];
                obj.Params.Tracula.in.fn = obj.Params.Eddy.out.fn{1} ;
                obj.Params.Tracula.in.dcmrirc = [obj.dependencies_dir  filesep 'TRACULA_DEPS' filesep  'dcmrirc.template' ];
                obj.Params.Tracula.in.FSDIR = obj.Params.FreeSurfer.dir;
                obj.Params.Tracula.in.bvec = obj.Params.Eddy.out.bvecs{1};
                obj.Params.Tracula.in.bval = obj.Params.Eddy.in.bvals{1};
                obj.Params.Tracula.in.nb0 = 5;
                %IMPLEMENTED BUT FULLY TESTED HENCE COMMENTED BELOW:
                obj.proc_tracula();
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %(INACTIVE) TRKLAND_HIPPOCING:
            for tohide=1:1
%             obj.Trkland.hippocing.in.hippo_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Left-Hippocampus.nii.gz');
%             obj.Trkland.hippocing.in.hippo_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Right-Hippocampus.nii.gz');
%             obj.Trkland.hippocing.in.postcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-posteriorcingulate.nii.gz');
%             obj.Trkland.hippocing.in.postcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-posteriorcingulate.nii.gz');
%             %Interpolation n (for cingulum):
%             obj.Trkland.hippocing.in.n_interp=33;
%             %DONE BUT ADDED TO TRKS_DATA (NEED TO QC FIRST HENCE COMMENTED
%             %BELOW: 
%             %obj.trkland_hippocing(); 

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %(INACTIVE) TRKLAND_CINGULUM:
            for tohide=1:1
%             obj.Trkland.cingulum.in.rostantcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-rostralanteriorcingulate.nii.gz');
%             obj.Trkland.cingulum.in.rostantcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-rostralanteriorcingulate.nii.gz');
%             obj.Trkland.cingulum.in.postcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-posteriorcingulate.nii.gz');
%             obj.Trkland.cingulum.in.postcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-posteriorcingulate.nii.gz');
%            
%             % BASE ON THESE VALUES (FROM EARLIER ADRC_PROCESSING W/O INTERP:
%             % RESOLUTION IN ADRC CONNECTOME DATA IS 1.8^3, here 2.0^3 so I
%             % decided to use the same number of interpolation points. 
%             % ninter_fx = 40;
%             % ninter_cingulum = 32;
%             % ninter_hippocing = 33;
%             obj.Trkland.cingulum.in.n_interp = 32;
%             %DONE BUT ADDED TO TRKS_DATA (NEED TO QC FIRST HENCE COMMENTED BELOW) 
            %trkland_cingulum(obj); obj.resave();
            end
        end
    
        function obj = setMyParams(obj)
            %software used:
            obj.dsistudio_version='GQI_DSISv053117';
            
            %root directoy where raw data lives:
            obj.object_dir= '/cluster/brutha/MATLAB_Scripts/Objects/Diffusion/'; %To add to path if needed
            obj.dcm_location='/cluster/sperling/HAB/Project1/DICOM_ARCHIVE/All_DICOMS/';
            obj.session_dir='/cluster/sperling/HAB/Project1/Sessions/';
            
            %template dependencies:
            obj.HABn272_meanFA='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/HABn272_MNI_Target/HABn272_meanFA.nii.gz';
            obj.HABn272_meanFA_skel_dst='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/HABn272_MNI_Target/HABn272_meanFA_skeleton_mask_dst.nii.gz';
            obj.ref_region='/usr/pubsw/packages/fsl/5.0.9/data/standard/LowerCingulum_1mm.nii.gz';
            
            %sh dependencies:
            obj.init_rotate_bvecs_sh='/cluster/brutha/MATLAB_Scripts/Objects/Diffusion/deps/SHELL/mod_fdt_rotate_bvecs.sh'; %For standarizing the bvecs after proc_dcm2nii
            obj.col2rows_sh='/cluster/bang/ADRC/Scripts/DEPENDENCIES/PREPROC_DEPS/drigo_col2rows.sh';
            obj.dependencies_dir='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/';
            %Freesurfer Dependencies:
            obj.init_FS = '/usr/local/freesurfer/stable6';
            obj.FS_location='/cluster/bang/HAB_Project1/FreeSurferv6.0';
            obj.FSL_dir = '/usr/pubsw/packages/fsl/5.0.9/';
            
            %DataCentral (if false, we won't upload)
            obj.dctl_flag = false;
            
            %skel_TOI dependencies
            obj.skeltoi_location='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/edJHU_MASK_ROIs/';
            obj.skeltoi_tois={ 'not_ATR_edJHU' 'not_CCG_edJHU' 'not_CHIPP_edJHU' 'not_CST_edJHU'  ...
                'not_Fma_edJHU' 'not_Fmi_edJHU' 'not_IFOF_edJHU' 'not_ILF_edJHU' ...
                'not_SLF_no_temp_edJHU'  'ATR_L_edJHU' 'ATR_R_edJHU' ...
                'CCG_L_edJHU' 'CCG_R_edJHU' 'CHIPP_L_edJHU' 'CHIPP_R_edJHU' 'CST_L_edJHU' ...
                'CST_R_edJHU' 'Fma_edJHU' 'Fmi_edJHU' 'IFOF_L_edJHU' 'IFOF_R_edJHU' ...
                'ILF_L_edJHU' 'ILF_R_edJHU' 'SLF_L_edJHU' 'SLF_L_edJHU_no_temporal' ...
                'SLF_R_edJHU' 'SLF_R_edJHU_no_temporal' 'allTracts_edJHU' 'Global_skel'  };
            
            %trkland dependencies:
            obj.fx_template_dir='/autofs/cluster/bang/ADRC/TEMPLATES/FX_1.8mm_orig/';
            obj.redo_history = false; %(always FALSE unles...) --> Allows us to redo the history of all processes withouth running any obj.BashCode.
            
            %Global parameters:
            obj.vox= [2 2 2 ];
            %obj.setDefaultParams; %from dwiMRI_Session class NOT USED!
            
        end
     
        function resave(obj)
            save(strrep([obj.objectHome filesep obj.sessionname '.mat'], [ filesep filesep ], filesep ),'obj');
        end
    end
    
    methods ( Access = protected )
        function obj = getDCM2nii(obj,torun)
            ii=1; %Expecting only on DWI sequence in the HAB1 project, so no loop
            
            %For proc_DCM2NII:
            obj.Params.DCM2NII.specific_vols=35;
            
            obj.Params.DCM2NII.scanlog = [ obj.session_location  filesep 'LogFiles' ...
                filesep 'scan.log' ] ;
            
            obj.Params.DCM2NII.newUnpack=false; %Flag that denotes whether the older and newer Unpack was used to dcm2nii the data.
            %Check early unpack dependency (unpacksdcmdir creates
            %obj.Params.DCM2NII.scanlog ) 
            if ~exist(obj.Params.DCM2NII.scanlog,'file')
                obj.Params.DCM2NII.scanlog = 'DNE'; %DOES NOT EXIST IN LATER NewUnpack.m
                obj.Params.DCM2NII.rawDiff = strrep(obj.root,'DWIs','Diffusion_HighRes');
                if exist(obj.Params.DCM2NII.rawDiff,'dir') ~=7
                    warning(['No scanlog found in:' obj.Params.DCM2NII.scanlog] ) ;
                    warning(['Or Diffusion_No scanlog found in:' obj.Params.DCM2NII.scanlog] ) ;
                else
                    obj.Params.DCM2NII.newUnpack=true;
                end
            end
            
            obj.Params.DCM2NII.seq_names='DIFFUSION_HighRes_30';
            %we modified grep err instead of tail -1 due to error in some subjects (e.g. 120419_4TT01420)
            if obj.Params.DCM2NII.newUnpack
                obj.Params.DCM2NII.out.rawOrig = strtrim(dir_wfp([obj.Params.DCM2NII.rawDiff '*.nii' ]));
                if iscell(obj.Params.DCM2NII.out.rawOrig)
                    obj.Params.DCM2NII.out.rawOrig=cell2char_rdp(obj.Params.DCM2NII.out.rawOrig);
                end
                exec_cmd=[ 'fslinfo ' obj.Params.DCM2NII.out.rawOrig   ' | grep ^dim4 | awk ''{ print $2 }'' ' ];
                [ ~ , temp_nvols ] = system(exec_cmd); 
                
                %Added the 'in' fields if it doesnt exist:
                if ~isfield(obj.Params.DCM2NII,'in')
                    obj.Params.DCM2NII.in = [] ; 
                end
                %Assigning number of volumes:
                if ischar(temp_nvols)
                    obj.Params.DCM2NII.in.nvols = str2num(temp_nvols);
                else %assuming double...
                    obj.Params.DCM2NII.in.nvols = temp_nvols;
                end
                
                %CHECK IF CORRECT NUMBER OF NVOLS:
                if obj.Params.DCM2NII.in.nvols ~=35
                    error(['This image does not contain 35 volume (for' obj.projectID '): ' obj.Params.DCM2NII.rawDiff  ])
                end
                
            else
                exec_cmd=[ 'cat ' obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names ' | grep " 35 " | grep err | awk ''{ print $7 }'' ' ];
                %we modified grep err instead of tail -1 due to error in some subjects (e.g. 120419_4TT01420)
                exec_cmd=[ 'cat ' obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names ' | grep " 35 " | grep err | awk ''{ print $8 }'' ' ];
                [ ~ , obj.Params.DCM2NII.in(ii).first_dcmfiles ] = system(exec_cmd);
            end
   
            
            obj.Params.DCM2NII.in.fsl2std_param = '-1 0 0 254 \n0 1 0 254 \n0 0 -1 0 \n0 0 0 1';
            obj.Params.DCM2NII.out(ii).location = [ obj.root 'Orig' filesep ];
            obj.Params.DCM2NII.out(ii).fn = [ obj.Params.DCM2NII.out(ii).location obj.Params.DCM2NII.seq_names '.nii.gz' ];
            
            
            if (torun) ; 
                obj.proc_dcm2nii() ; 
            end
        end            
    end
end