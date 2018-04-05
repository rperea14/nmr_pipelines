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
        %software used:
        dsistudio_version='GQI_DSISv053117'
        
        %root directoy where raw data lives:
        object_dir= '/cluster/brutha/MATLAB_Scripts/Objects/Diffusion/'; %To add to path if needed
        root_location='/cluster/sperling/HAB/Project1/Sessions/';
        dcm_location='/cluster/sperling/HAB/Project1/DICOM_ARCHIVE/All_DICOMS/';
        session_location='/cluster/sperling/HAB/Project1/Sessions/';
        
        %template dependencies:
        HABn272_meanFA='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/HABn272_MNI_Target/HABn272_meanFA.nii.gz';
        HABn272_meanFA_skel_dst='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/HABn272_MNI_Target/HABn272_meanFA_skeleton_mask_dst.nii.gz';
        ref_region='/usr/pubsw/packages/fsl/5.0.9/data/standard/LowerCingulum_1mm.nii.gz';
        
        %sh dependencies:
        init_rotate_bvecs_sh='/cluster/brutha/MATLAB_Scripts/Objects/Diffusion/deps/SHELL/mod_fdt_rotate_bvecs.sh'; %For standarizing the bvecs after proc_dcm2nii
        col2rows_sh='/cluster/bang/ADRC/Scripts/DEPENDENCIES/PREPROC_DEPS/drigo_col2rows.sh';
        dependencies_dir='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/';
        %Freesurfer Dependencies:
        init_FS = '/usr/local/freesurfer/stable6';
        FS_location='/cluster/bang/HAB_Project1/FreeSurferv6.0';
        %DataCentral (if false, we won't upload)
        dctl_flag = false;
        
        %skel_TOI dependencies
        skeltoi_location='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/edJHU_MASK_ROIs/';
        skeltoi_tois={ 'not_ATR_edJHU' 'not_CCG_edJHU' 'not_CHIPP_edJHU' 'not_CST_edJHU'  ...
            'not_Fma_edJHU' 'not_Fmi_edJHU' 'not_IFOF_edJHU' 'not_ILF_edJHU' ...
            'not_SLF_no_temp_edJHU'  'ATR_L_edJHU' 'ATR_R_edJHU' ...
            'CCG_L_edJHU' 'CCG_R_edJHU' 'CHIPP_L_edJHU' 'CHIPP_R_edJHU' 'CST_L_edJHU' ...
            'CST_R_edJHU' 'Fma_edJHU' 'Fmi_edJHU' 'IFOF_L_edJHU' 'IFOF_R_edJHU' ...
            'ILF_L_edJHU' 'ILF_R_edJHU' 'SLF_L_edJHU' 'SLF_L_edJHU_no_temporal' ...
            'SLF_R_edJHU' 'SLF_R_edJHU_no_temporal' 'allTracts_edJHU' 'Global_skel'  };
        
        %trkland dependencies:
        fx_template_dir= '/space/public_html/rdp20/fornix_ROA/FX_1.8mm_orig/';
        redo_history = false; %(always FALSE unles...) --> Allows us to redo the history of all processes withouth running any obj.BashCode.
    end
    
    methods
        function obj = dwi_HAB(sessionname,opt)
            %For compiler code:
            if ~isdeployed()
                addpath(genpath(obj.object_dir));
            end
            %%%  If opt is passed, then the root Sessions folder will be
            %%%  replaced with this argument.
            if nargin>1
                if strcmp(opt,'DataCentral')
                    obj.dctl_flag = true ; 
                else
                    obj.root = opt;
                end
            else
                obj.dctl_flag = false ; 
            end
            
            %Check if you are in the right project (if not, probably
            %obj.session_location doesn't exist:
            if ~exist([obj.session_location sessionname],'dir')
                error(['Sesion directory: ' obj.session_location ' doesn''t exist. Are you sure you are in the right Project ID?']);
            end
            
            %Initialize root variables:
            obj.sessionname = sessionname;
            obj.root = [obj.root_location sessionname '/DWIs/'];
            obj.dcm_location= [ obj.dcm_location sessionname filesep ];
            obj.session_location= [ obj.session_location sessionname filesep ];
            
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
            
            %Assign Project ID:
            obj.projectID='HAB';
            
            %Check if *.nii.gz files exist, if not get them from DCM2nii:
            obj.rawfiles = dir_wfp([obj.root 'Orig/*.nii.gz' ] );
            if isempty(obj.rawfiles)
                RunFlag=true;
                obj.getDCM2nii(RunFlag);
            end
            
            if nargin>1
                if strcmp(opt,'DataCentral')
                    obj.dctl_flag = true ;
                    obj.proc_get_eddymotion();
                    obj.UploadData_DWI();
                    return
                else
                    if ~strcmpi(oldroot,newroot)
                        obj = replaceObjText(obj,{oldroot},{newroot});
                        obj.resave;
                    end
                end
            else
                obj.dctl_flag = false ;
            end
            %Continue with CommonPreProc
            obj.CommonPreProc();
            %Continue with CommonPostProc
            obj.CommonPostProc();
            
        end
        
        function obj = setMyParams(obj)
            %Global parameters:
            obj.vox= [2 2 2 ];
            obj.setDefaultParams; %from dwiMRI_Session class
        end
        
        function obj = CommonPreProc(obj)
            obj.dosave = true ; %To record process in MAT file
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Get DCM2NII File location:
            %
            RunFlag=false; %here, we only allocate variable, not run proc_dcm2nii
            obj.getDCM2nii(RunFlag);
            
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For FreeSurfer recon-all
            [ tmpa, tmpb ] = system('whoami ');
            %[ tmpa, tmpb ] = system('echo $0');
            %[~,  tmpshell , ~] = fileparts(tmpb);
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
            
            %A possible error is the naming convention when only a T1 was
            %used!!
              
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
            
            obj.proc_T1toDWI();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if obj.dctl_flag == true
                if ~isdeployed() %Not compiled code, so run this!
                    %Uploading Skel Data into DataCentral:
                    UploadData_DWI(obj)
                    
                    if obj.dosave
                        save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
                    end
                end
            end
        end
        
        
        function obj = CommonPostProc(obj) 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Creating Fornix TRKLAND 
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
                obj.trkland_fx(); 
               
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %TRKLAND_HIPPOCING:
            for tohide=1:1
            obj.Trkland.hippocing.in.hippo_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Left-Hippocampus.nii.gz');
            obj.Trkland.hippocing.in.hippo_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Right-Hippocampus.nii.gz');
            obj.Trkland.hippocing.in.postcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-posteriorcingulate.nii.gz');
            obj.Trkland.hippocing.in.postcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-posteriorcingulate.nii.gz');
            %Interpolation n (for cingulum):
            obj.Trkland.hippocing.in.n_interp=33;
            
            %DONE BUT ADDED TO TRKS_DATA (NEED TO QC FIRST HENCE COMMENTED
            %BELOW: 
            %obj.trkland_hippocing(); 
            %obj.resave();
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %TRKLAND_CINGULUM:
            for tohide=1:1
            obj.Trkland.cingulum.in.rostantcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-rostralanteriorcingulate.nii.gz');
            obj.Trkland.cingulum.in.rostantcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-rostralanteriorcingulate.nii.gz');
            obj.Trkland.cingulum.in.postcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-posteriorcingulate.nii.gz');
            obj.Trkland.cingulum.in.postcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-posteriorcingulate.nii.gz');
           
            % BASE ON THESE VALUES (FROM EARLIER ADRC_PROCESSING W/O INTERP:
            % RESOLUTION IN ADRC CONNECTOME DATA IS 1.8^3, here 2.0^3 so I
            % decided to use the same number of interpolation points. 
            % ninter_fx = 40;
            % ninter_cingulum = 32;
            % ninter_hippocing = 33;
            obj.Trkland.cingulum.in.n_interp = 32;
            %DONE BUT ADDED TO TRKS_DATA (NEED TO QC FIRST HENCE COMMENTED BELOW) 
            %trkland_cingulum(obj); obj.resave();
            end
            
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %TRACULA (and implicit functionality of bedpostx):
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
        function resave(obj)
            save(strrep([obj.objectHome filesep obj.sessionname '.mat'], [ filesep filesep ], filesep ),'obj');
        end
    end
    
    methods ( Access = protected )
        function obj = getDCM2nii(obj,torun)
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
                exec_cmd=[ 'TEMP_AA=$( ls ' obj.Params.DCM2NII.rawDiff '*.nii ) ; fslinfo $TEMP_AA | grep ^dim4 | awk ''{ print $2 }'' ' ];
            else
                exec_cmd=[ 'cat ' obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names ' | grep " 35 " | grep err | awk ''{ print $7 }'' ' ];
                %we modified grep err instead of tail -1 due to error in some subjects (e.g. 120419_4TT01420)
                exec_cmd=[ 'cat ' obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names ' | grep " 35 " | grep err | awk ''{ print $8 }'' ' ];
                [ ~ , obj.Params.DCM2NII.in(ii).first_dcmfiles ] = system(exec_cmd);
            end
            [ ~ , obj.Params.DCM2NII.in.nvols ] = system(exec_cmd);
            
            obj.Params.DCM2NII.in.nvols=str2num(strtrim(obj.Params.DCM2NII.in.nvols));
            %CHECK IF CORRECT NUMBER OF NVOLS:
            if obj.Params.DCM2NII.in.nvols ~=35 
                error(['Raw image: '])
            end
            
            ii=1; %Expecting only on DWI sequence in the HAB1 project, so no loop
            obj.Params.DCM2NII.in.fsl2std_param = '-1 0 0 254 \n0 1 0 254 \n0 0 -1 0 \n0 0 0 1';
            obj.Params.DCM2NII.out(ii).location = [ obj.root 'Orig' filesep ];
            obj.Params.DCM2NII.out(ii).fn = [ obj.Params.DCM2NII.out(ii).location obj.Params.DCM2NII.seq_names '.nii.gz' ];
            
            
            if (torun) ; 
                obj.proc_dcm2nii() ; 
            end
        end            
    end
    
    methods
        %UNDER DEVELOPMENT
           %UNDER DEVELOPMENT
           function obj = post_tracx_by_txt(obj,masktxt_fname,masktxt_dir,replace_masktxt_info)
               %%%%%% CODE FOR DEALING WITH DIFFERNT MASK FOR PROB TRACTOGRAPHY%
               for tohide=1:1
                   if nargin<2 || isempty(masktxt_fname)
                       masktxt_dir=[obj.dependencies_dir 'fMRI_masks' filesep 'mask_txt' filesep ];
                       masktxt_fname = 'default_mask'; %obj.dep_dir / fMRI_masks/maskt_txt/default_mask.txt
                   end
                   if nargin < 3 || isempty(masktxt_dir)
                       masktxt_dir=[obj.dependencies_dir 'fMRI_masks' filesep 'mask_txt' filesep ];
                   end
                   
                   %Verify you don't want to replace the processing of the
                   %specific mask:
                   if nargin <4
                       replace_masktxt_info=false;
                   end
                   
                   %Initialize the Params structure if not already:
                   if ~isfield(obj.Params,'tracxBYmask')
                       obj.Params.tracxBYmask=[];
                   end
                   
                   %Check if this is the first txt_filename that is inputted
                   %*This will allow us to run different instances of tractxBYmask
                   %without losing information about what was ran:
                   [~ , tmp_txtfname, ~ ]  = fileparts(masktxt_fname);
                   if ~isfield(obj.Params.tracxBYmask,'list_txt_fnames')
                       obj.Params.tracxBYmask.list_txt_fnames{1} = {tmp_txtfname};
                   else
                       %Double check that the txt_fname hasn't been used or else
                       %it will be replaced
                       flag_list_exist=0;
                       for ijk=1:numel(obj.Params.tracxBYmask.list_txt_fnames)
                           if strcmp(tmp_txtfname,obj.Params.tracxBYmask.list_txt_fnames{ijk})
                               if replace_masktxt_info ~= 1
                                   warning(['The name for the mask_txt fname: ' tmp_txtfname ' has been used (and probably already ran, not necessarily successfully).'])
                                   display(['Either 1) change the masktxt_fname (1st argument in obj.post_tracx_by_txt() ) filename or' ])
                                   display('       2) check the obj.Params.tracxBYmask.(masktxt_fname) structure or')
                                   display('       3) re-run the obj.post_tracx_by_txt with a 3rd argument as True (*info data will be replaced)')
                                   display('Returning...')
                                   return
                               end
                               %This flag tells me that the same name exists!
                               flag_list_exist=1;
                           end
                       end
                       %Add the newer txt_fname is the list of tracx parameters,
                       %it replace_masktxt_info is true then no need as the values
                       %will be explicitly be replace
                       if flag_list_exist == 0
                           obj.Params.tracxBYmask.list_txt_fnames{end+1} = {tmp_txtfname};
                       else
                           warning(['All the Params info for obj.Params.tracxBTmask.' tmp_txtfname ' will be replaced...'])
                       end
                   end
               end
               %%%%%% END DEALING WITH VARIABLE NUMBER OF MASK FOR TRACX%%%%%%
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
               %%%%%%%%%%%%%%%%%% CHECKING ARGUMENTS INPUTTED%%%%%%%%%%%%%%%%%
               obj.Params.tracxBYmask.allmasks.(tmp_txtfname).in.dir = masktxt_dir;
               for tohide=1:1
                   %Check to see whether a mask directory locations is input
                   %Check to see what mask_txt file is inputted
                   if nargin <2
                       tmp_txt_fullpath = [ obj.Params.tracxBYmask.allmasks.(tmp_txtfname).in.dir 'try_masks.txt' ] ;
                       display(['DENOTING: ' tmp_txt_fullpath ' as the default mask filename to process (no arguments passed).'])
                       pause(1)
                   else
                       tmp_txt_fullpath =[ obj.Params.tracxBYmask.allmasks.(tmp_txtfname).in.dir masktxt_fname '.txt' ] ;
                       display(['DEFAULT MASK FILENAME: ' tmp_txt_fullpath])
                   end
                   
                   %
               end
               %%%%%%%%%%%%%%%END CHECKING ARGUMENTS INPUTTED%%%%%%%%%%%%%%%%%
               
               
               %%%%%%%%%%%%%%%%%% VARIABLE INITIALIZATION%%%%%%%%%%%%%%%%%%%%%
               %VARIABLE INITIALIZATION:
               %TRACULA related:
               obj.Params.tracxBYmask.tracula.bedp_dir = fileparts(obj.Params.Tracula.out.bedp_check);
               obj.Params.tracxBYmask.T1_tmp = [fileparts(which('spm')) '/canonical/single_subj_T1.nii'];;
               obj.Params.tracxBYmask.tracula.b0 = [obj.Params.tracxBYmask.tracula.bedp_dir ...
                   filesep '..' filesep 'dmri' filesep 'lowb.nii.gz' ];
               
               
               %TXT file related:
               obj.Params.tracxBYmask.allmasks.(tmp_txtfname).in.txt_fname = tmp_txt_fullpath;
               
               
               obj.Params.tracxBYmask.allmasks.(tmp_txtfname).in.movefiles = [ '.' filesep 'post_tracx' filesep 'all_masks' filesep  tmp_txtfname ];
               %obj.Params.tracxBYmask.allmasks.(tmp_txtfname).in.movefiles = ...
               %      ['..' filesep '..' filesep '..' filesep 'post_tracx' filesep 'all_masks' filesep  tmp_txtfname ];
               
               obj.Params.tracxBYmask.allmasks.(tmp_txtfname).probtracx2_args = ...
                   ' -l --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd  ' ;
               
               proc_tracxBYmask(obj,tmp_txtfname); %obj.resave()
               %%%%%%%%%%%%%%% END VARIABLE INITIALIZATION%%%%%%%%%%%%%%%%%%%%
               
               
               %%%%%%%%%%%%%%%% IMPLEMENTATION STARTS HERE %%%%%%%%%%%%%%%%%%%
               %%%%%%%%%%%%%%%%%% END OF IMPLEMENTATION  %%%%%%%%%%%%%%%%%%%%%
               
           end
    end
end