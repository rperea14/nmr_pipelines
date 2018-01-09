classdef dwi_HAB_beta < dwiMRI_Session_beta
    %%  classdef dwi_ADRC < dwiMRI_Session_beta
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
        root_location='/cluster/sperling/HAB/Project1/DWIs_30b700/Sessions/';
        dcm_location='/cluster/sperling/HAB/Project1/DICOM_ARCHIVE/All_DICOMS/';
        session_location='/cluster/sperling/HAB/Project1/Sessions/';
        
        %template dependencies:
        HABn272_meanFA='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/HABn272_MNI_Target/HABn272_meanFA.nii.gz';
        HABn272_meanFA_skel_dst='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/HABn272_MNI_Target/HABn272_meanFA_skeleton_mask_dst.nii.gz';
        ref_region='/usr/pubsw/packages/fsl/5.0.9/data/standard/LowerCingulum_1mm.nii.gz';
       
        %sh dependencies:
        init_rotate_bvecs_sh='/cluster/sperling/HAB/Project1/Scripts/DWIs/mod_fdt_rotate_bvecs.sh'; %For standarizing the bvecs after proc_dcm2nii
        dependencies_dir='/cluster/hab/HAB/Project1/DWIs_30b700/DEPENDENCIES/';
        %Freesurfer Dependencies:
        init_FS = '/usr/local/freesurfer/stable6';

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
        
    end
    
    methods
        function obj = dwi_HAB_beta(sessionname,opt)
            
            %For compiler code:
            if ~isdeployed()
                addpath(genpath('/autofs/space/kant_004/users/rdp20/scripts/matlab'));
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
                obj.setMyParams;
            end
            
            %Check if *.nii.gz files exist, if not get them from DCM2nii:
            obj.rawfiles = dir_wfp([obj.root 'Orig/*.nii.gz' ] );
            if isempty(obj.rawfiles)
                RunFlag=true;
                obj.getDCM2nii(RunFlag);
            end
            
            if nargin>1
                if strcmp(opt,'DataCentral')
                    obj.dctl_flag = true ;
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
            obj.setDefaultParams; %from dwiMRI_Session_beta class
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
            %To calculate mean motion based on edddy output:
            obj.Params.EddyMotion.in.movefiles = ['..' filesep '03_EddyMotion']; 
            obj.Params.EddyMotion.in.fn_eddy = obj.Params.Eddy.out.fn ;
            
            obj.proc_get_eddymotion();
                        
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For B0mean:
            obj.Params.B0mean.in.movefiles = ['..' filesep '03_B0mean'];
            obj.Params.B0mean.in.fn=obj.Params.Eddy.out.fn;
            obj.Params.B0mean.in.b0_nvols=5;
            
            obj.proc_meanb0();
            
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
            obj.Params.GQI.in.bvecs = obj.Params.Eddy.out.bvecs;
            obj.Params.GQI.in.bvals = obj.Params.Eddy.in.bvals;
            obj.Params.GQI.in.mask = obj.Params.MaskAfterEddy.out.finalmask;
            obj.Params.GQI.in.prefix = obj.dsistudio_version ; %Double check this so you prefix the version of DSISTUDIO!
            obj.Params.GQI.out.export = 'gfa,nqa0,nqa1';
            
            obj.proc_gqi();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For AntsReg:
            obj.Params.AntsReg.in.movefiles = ['..' filesep '06_Ants_CoReg' ];
            obj.Params.AntsReg.in.fn = obj.Params.Dtifit.out.FA ;
            obj.Params.AntsReg.in.ref = obj.HABn272_meanFA;
            obj.Params.AntsReg.in.threads = '4' ;
            obj.Params.AntsReg.in.prefix = 'Antsv201_2_HABn272_' ;
            
            obj.proc_antsreg();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For Skeletonize:
            obj.Params.Skeletonize.in.movefiles = ['..' filesep '07_Skeletonize' ];
            obj.Params.Skeletonize.in.fn = obj.Params.AntsReg.out.fn ;
            obj.Params.Skeletonize.in.meanFA = obj.HABn272_meanFA;
            obj.Params.Skeletonize.in.skel_dst = obj.HABn272_meanFA_skel_dst;
            obj.Params.Skeletonize.in.thr = '0.3';
            obj.Params.Skeletonize.in.ref_region = obj.ref_region;
            obj.Params.Skeletonize.in.prefix = [ 'FSLv509_skelHABn272' ] ; %check fsl versio you are calling!
            
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
            if obj.dctl_flag == true
                if ~isdeployed() %Not compiled code, so run this!
                    %Uploading Skel Data into DataCentral:
                    UploadData_DWI(obj)
                    
                    if obj.dosave
                        save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
                    end
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For FreeSurfer recon-all
            [ tmpa, tmpb ] = system('whoami ');
            %[ tmpa, tmpb ] = system('echo $0');
            %[~,  tmpshell , ~] = fileparts(tmpb);
            obj.Params.FreeSurfer.shell = strtrim(tmpb); %strtrim(tmpshell);
          
            obj.Params.FreeSurfer.dir = [ filesep 'eris' filesep 'bang' filesep ...
            'HAB_Project1' filesep 'FreeSurferv6.0' filesep ] ;
            obj.Params.FreeSurfer.init_location = obj.init_FS; 
            %Retrieving a T1 scan:
            [sys_error, obj.Params.FreeSurfer.in.T1 ] = system(['ls ' obj.session_location 'MPRAGE' filesep '*.mgz | head -1' ]);
            if sys_error ~= 0 %No problem, we get the T1 the continue...
                fprintf(['\nError when finding the T1:'  obj.Params.FreeSurfer.in.T1  '\n'])
            end
            
            %Retrieving a T2 scan:
            [sys_error, obj.Params.FreeSurfer.in.T2 ] = system(['ls ' obj.session_location 'other' filesep '*T2* | head -1' ]);
            if sys_error ~= 0 %No problem, we get the T1 the continue...
                fprintf(['\nNo T2 found:'  obj.Params.FreeSurfer.in.T2  '\n'])
                obj.Params.FreeSurfer.in.T2exist=false;
            else
                obj.Params.FreeSurfer.in.T2exist=true;
            end
            
            obj.Params.FreeSurfer.out.aparcaseg = [ obj.Params.FreeSurfer.dir ...
                filesep obj.sessionname filesep 'mri' filesep 'aparc+aseg.mgz' ] ;
            
            %obj.proc_getFreeSurfer();
             
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %FS2dwi (move aparc and aparc2009 segmes to dwi space):
            obj.Params.FS2dwi.in.movefiles = ['..' filesep '05_FS2dwi' ];
            obj.Params.FS2dwi.in.b0 = obj.Params.B0mean.out.fn ; 
            obj.Params.FS2dwi.in.aparcaseg = obj.Params.FreeSurfer.out.aparcaseg ; 
            
            obj.Params.FS2dwi.in.tmpfile_aparcaseg = [ obj.dependencies_dir 'FS_aparc.txt' ] ; 
            obj.Params.FS2dwi.in.tmpfile_aparcaseg2009 = [ obj.dependencies_dir 'FS_aparc2009.txt' ] ; 
            obj.Params.FS2dwi.in.tmpfile_hippo_bil = [ obj.dependencies_dir 'FS_hippolabels_bil.txt' ] ;
            
            
            
            obj.Params.FS2dwi.in.aparcaseg2009 = ...
                strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','aparc.a2009s+aseg')); 
            
            %A possible error is the naming convention when only a T1 was
            %used!!
            obj.Params.FS2dwi.in.hippofield_left = ...
                strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','lh.hippoSfLabels-T1-T2.v10.FSvoxelSpace')); 
            obj.Params.FS2dwi.in.hippofield_right = ...
                strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','rh.hippoSfLabels-T1-T2.v10.FSvoxelSpace')); 
            
            obj.proc_FS2dwi();
            
        end
        
        function obj = CommonPostProc(obj) 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %TRACULA (and implicit functionality of bedpostx):
            obj.Params.Tracula.in.movefiles = ['..' filesep 'post_TRACULA' ];
            obj.Params.Tracula.in.fn = obj.Params.Eddy.out.fn{1} ; 
            obj.Params.Tracula.in.dcmrirc = [obj.dependencies_dir 'dcmrirc.template' ];
            obj.Params.Tracula.in.FSDIR = obj.Params.FreeSurfer.dir;
            obj.Params.Tracula.in.bvec = obj.Params.Eddy.out.bvecs{1};  
            obj.Params.Tracula.in.bval = obj.Params.Eddy.in.bvals{1};
            obj.Params.Tracula.in.nb0 = 5;
            obj.Params.Tracula.in.prefix = 'hab';
            
        %    obj.proc_tracula();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Creating Fornix TRKLAND 
            %TRKLAND
            obj.Trkland.root = [ obj.root  'post_TRKLAND' filesep ];
            %FX_TRKLAND:
            for tohide=1:1
                %b0 params:
                obj.Trkland.fx.in.b0 = obj.Params.B0mean.out.fn{end};
                %Template parameters:
                obj.Trkland.fx.tmp.b0 = [ obj.fx_template_dir '141210_8CS00178_b0.nii.gz' ] ;
                obj.Trkland.fx.tmp.roa_solid_bil =[ obj.fx_template_dir 'TMP_178_bil_fx_dil11.nii.gz' ] ;
                obj.Trkland.fx.tmp.roa_solid_lh = [ obj.fx_template_dir 'TMP_178_lh_fx_dil11.nii.gz' ] ;
                obj.Trkland.fx.tmp.roa_solid_rh = [ obj.fx_template_dir 'TMP_178_rh_fx_dil11.nii.gz' ] ;
                obj.Trkland.fx.tmp.roi_bil = [ obj.fx_template_dir 'TMP_178_bil_fx_dil.nii.gz' ] ;
                obj.Trkland.fx.tmp.roi_lh = [ obj.fx_template_dir 'TMP_178_lh_fx_dil.nii.gz' ] ;
                obj.Trkland.fx.tmp.roi_rh = [ obj.fx_template_dir 'TMP_178_rh_fx_dil.nii.gz' ] ;
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
                
                obj.trkland_fx();
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %TRKLAND_HIPPOCING:
            for tohide=1:1
            obj.Trkland.hippocing.in.hippo_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Left-Hippocampus.nii.gz');
            obj.Trkland.hippocing.in.hippo_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Right-Hippocampus.nii.gz');
            obj.Trkland.hippocing.in.postcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-posteriorcingulate.nii.gz');
            obj.Trkland.hippocing.in.postcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-posteriorcingulate.nii.gz');
            
            obj.trkland_hippocing();
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %TRKLAND_CINGULUM:
            for tohide=1:1
            obj.Trkland.cingulum.in.rostantcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-rostralanteriorcingulate.nii.gz');
            obj.Trkland.cingulum.in.rostantcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-rostralanteriorcingulate.nii.gz');
            obj.Trkland.cingulum.in.postcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-posteriorcingulate.nii.gz');
            obj.Trkland.cingulum.in.postcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-posteriorcingulate.nii.gz');
            
            trkland_cingulum(obj)
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %AFQ:
             obj.Params.AFQ.in.movefiles = ['..' filesep 'post_AFQ' ];
            obj.Params.AFQ.in.dwi = obj.Params.Eddy.out.fn{1};
            
            obj.Params.AFQ.in.T1 = strtrim(obj.Params.FreeSurfer.in.T1);
            obj.Params.AFQ.in.tmp_bvecs = obj.Params.Eddy.out.bvecs{1};
            obj.Params.AFQ.in.tmp_bvals = obj.Params.Eddy.in.bvals{1};
            
            
            %Parameter specification for pre-processedAFQ
            obj.Params.AFQ.dwParams.bvalue = []; 
            obj.Params.AFQ.dwParams.gradDirsCode = [];
            obj.Params.AFQ.dwParams.clobber = 1; % 1 -> files overwritten, 0 -> ask for it, -1 -> existing files will be used
            obj.Params.AFQ.dwParams.dt6BaseName = ''; %if empty dti<nDirs>trilin
            obj.Params.AFQ.dwParams.flipLrApFlag =false ; 
            obj.Params.AFQ.dwParams.numBootStrapSamples  = 500; 
            obj.Params.AFQ.dwParams.fitMethod = 'ls';
            obj.Params.AFQ.dwParams.nStep = 50;
            obj.Params.AFQ.dwParams.eddyCorrect = 0; % 1 --> eddy & motion, 0 --> motion -1 --> nothing
            obj.Params.AFQ.dwParams.excludeVols = [];
            obj.Params.AFQ.dwParams.bsplineInterpFlag = false; %false = trilinear, true = bspline
            obj.Params.AFQ.dwParams.phaseEncodeDir = 2;
            obj.Params.AFQ.dwParams.dwOutMm = [1.8 1.8 1.8];
            obj.Params.AFQ.dwParams.rotateBvecsWithRx       = false;
            obj.Params.AFQ.dwParams.rotateBvecsWithCanXform = false;
            obj.Params.AFQ.dwParams.bvecsFile =  '';%initialize inside the superclass obj.Params.CoRegMultiple.out.combined_bvals;
            obj.Params.AFQ.dwParams.bvalsFile = ''; %initialize inside the superclass 
            obj.Params.AFQ.dwParams.noiseCalcMethod = 'b0';
            obj.Params.AFQ.dwParams.outDir = ''; %this will initialize within the superclass method
%             
%             
            obj.proc_AFQ();
            
            
            
        end
        
        function resave(obj)
            save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
        end
    end
    
    methods ( Access = protected )
        function obj = getDCM2nii(obj,torun)
            %For proc_DCM2NII:
            obj.Params.DCM2NII.specific_vols=35;
            obj.Params.DCM2NII.scanlog = [ obj.session_location  filesep 'LogFiles' ...
                filesep 'scan.log' ] ;
            
            if ~exist(obj.Params.DCM2NII.scanlog,'file')
                error(['No scanlog found in:' obj.Params.DCM2NII.scanlog '. Exiting...']);
            end
            
            obj.Params.DCM2NII.seq_names='DIFFUSION_HighRes_30';
            try
                %we modified grep err instead of tail -1 due to error in some subjects (e.g. 120419_4TT01420)
                exec_cmd=[ 'cat ' obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names ' | grep " 35 " | grep err | awk ''{ print $7 }'' ' ];
                [ ~ , obj.Params.DCM2NII.in.nvols ] = system(exec_cmd);
            catch
                errormsg=['DCM2NII: No 35 vols. when reading scanlog located in: ' ...
                    obj.Params.DCM2NII.scanlog '\n' ];
                obj.UpdateErrors(errormsg);
            end
            obj.Params.DCM2NII.in.fsl2std_param = '-1 0 0 254 \n0 1 0 254 \n0 0 -1 0 \n0 0 0 1';
            for ii=1:1 %For object compatiblity with ADRC that contains 3 DWIs sequences
                obj.Params.DCM2NII.in(ii).nvols=str2double(obj.Params.DCM2NII.in.nvols);
                %we modified grep err instead of tail -1 due to error in some subjects (e.g. 120419_4TT01420)
                exec_cmd=[ 'cat ' obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names ' | grep " 35 " | grep err | awk ''{ print $8 }'' ' ];
                [ ~ , obj.Params.DCM2NII.in(ii).first_dcmfiles ] = system(exec_cmd);
                
                obj.Params.DCM2NII.out(ii).location = [ obj.root 'Orig' filesep ];
                obj.Params.DCM2NII.out(ii).fn = [ obj.Params.DCM2NII.out(ii).location obj.Params.DCM2NII.seq_names '.nii.gz' ];
                
            end
            if (torun) ; obj.proc_dcm2nii ; end
        end
    end
end