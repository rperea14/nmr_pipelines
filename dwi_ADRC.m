classdef dwi_ADRC < dwiMRI_Session
    %%  classdef dwi_ADRC < dwiMRI_Session
    %%  This class is a subclass of its parent class dwi_MRI_Session.m
    %%  (where it will inherent other methods).
    %%  Created by:
    %%              Rodrigo D. Perea grandrigo@gmail.com
    %%
    
    properties
        %root directoy where raw data lives:
        session_location='/eris/bang/ADRC/Sessions/';
        dcm_location = '/eris/bang/ADRC/DICOM_Archive/';
        gradfile='/autofs/space/kant_004/users/ConnectomeScanner/Scripts/adrc_diff_prep/bash/gradient_nonlin_unwarp/gradient_coil_files/coeff_AS302.grad';
        dependencies_dir='/eris/bang/ADRC/Scripts/DEPENDENCIES/';
        
        %sh dependencies:
        sh_gradfile=[ '/eris/bang/ADRC/Scripts/DEPENDENCIES/GradNonLin_Correc/run_mris_gradient_nonlin__unwarp_volume__batchmode_ADRC_v3.sh ' ...
            '/usr/pubsw/common/matlab/8.5 '];
        b0MoCo_rotate_bvecs_sh='/eris/bang/ADRC/Scripts/DEPENDENCIES/PREPROC_DEPS/rotate_bvecs.sh'; %For rotating the bvecs after proc_b0MoCo
        init_rotate_bvecs_sh='/eris/bang/ADRC/Scripts/DEPENDENCIES/PREPROC_DEPS/mod_fdt_rotate_bvecs.sh'; %For standarizing the bvecs after proc_dcm2nii
        col2rows_sh='/eris/bang/ADRC/Scripts/DEPENDENCIES/PREPROC_DEPS/drigo_col2rows.sh';
        
 
        %FreeSurfer Dependencies
        FS_location='/eris/bang/ADRC/FreeSurferv6.0/';
        init_FS = '/usr/local/freesurfer/stable6';
        
        %trkland dependencies:
        fx_template_dir= '/space/public_html/rdp20/fornix_ROA/FX_1.8mm_orig/';
        
        %frois dependencies
        FROIS_dir = '/eris/bang/ADRC/TEMPLATES/FROIS/'
        
        %Related to T1 (maybe optional?)
        bb = [-78 -112 -70; 78 76 90];
        T1_vox = [ 1 1 1 ] ; 
    end
    
    
    methods
        function obj = dwi_ADRC(sessionname,opt)
             %For compiler code:
            if ~isdeployed()
                addpath(genpath('/autofs/space/kant_004/users/rdp20/scripts/matlab'));
            end
            
            %%%  If opt is passed, then the root Sessions folder will be
            %%%  replaced with this argument.
            if nargin>1
                obj.root = opt;
            end
            obj.sessionname = sessionname;
            obj.root = [obj.session_location sessionname '/DWIs/'];
            obj.dcm_location = [ obj.dcm_location sessionname filesep ];
            obj.session_location= [ obj.session_location sessionname filesep ] ;
          
            %If the folder /DWIs/ does not exist, then create it!
            if exist(obj.root,'dir')==0
                obj.make_root();
            end
            obj.objectHome = obj.root ; 
            %!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %!!!!!!!
      
            %Check to see if a *.mat file exists.
            if exist([obj.objectHome filesep sessionname '.mat'],'file')>0
                load([obj.objectHome filesep sessionname '.mat']);
                oldroot = obj.root;
                obj.wasLoaded = true;
            else
              %  obj.setMyParams;
              donothing=1;
            end
            
            %Init project ID
            obj.projectID='ADRC';
            
            if isempty(obj.T1)
                  [~ , obj.T1 ]  = system([ 'ls ' obj.session_location 'T1' filesep '*.nii | head -1' ]);
            end
            
            %Check if *.nii.gz files exist, if not get them from DCM2nii:
            obj.rawfiles = dir_wfp([obj.root 'Orig/*.nii.gz' ] );
            if isempty(obj.rawfiles) || numel(obj.rawfiles) ~= 4 % 4 DWIs sequence acquired here
                RunFlag = true ;
                obj.getDCM2nii(RunFlag);
            end
        
            %Assign whether a single process is given or a particular one
            if nargin >1 
                obj.exec_onecmd = opt;
                obj.ParticularProc();
            else
                obj.exec_onecmd = '' ;
                %Start the CommonPreProc:
                obj.CommonPreProc();
                %Start the CommonPostProc (commenting it? Maybe, maybe not?):
                obj.CommonPostProc();
                fprintf(['\n>>*If first time running, please run obj.CommonPostProc() after this! \n<<\n'])
            end
        end
        
        function obj=setMyParams(obj)
            %%%%%%%%%%%%
            %Global parameters:
%             obj.vox = [1.8 1.8 1.8];
%             obj.setDefaultParams; %this will call the method in the superclass dwiMRI_Session.m 
%             obj.rawfiles = dir_wfp([obj.root 'Orig' filesep '*.nii.gz' ] );
%             
        end
        
         
        function resave(obj)
            save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
        end
        
        function obj = CommonPreProc(obj)
            obj.dosave = true ; %To record process in MAT file
            
             if isempty(obj.rawfiles)
                 obj.rawfiles = dir_wfp([obj.root 'Orig_' filesep '*.nii' ] );
             end
             
            
            %%%%%%%%%%%%
            %01_DropVols
            %For proc_dropvols
            obj.Params.DropVols.in.tmin='1';
            obj.Params.DropVols.in.tsize='67';
            obj.Params.DropVols.in.prefix='dv_';
            obj.Params.DropVols.in.movefiles=['..' filesep '01_DropVols' filesep ];
            obj.Params.DropVols.in.fn=obj.rawfiles;
            %Bvecs and bvals will be created from XX.in.fn and XX.out.fn
            obj.Params.DropVols.out.fn=dir_wfp([obj.root, '01_DropVols', filesep, '*.nii.gz']);
            
            obj.proc_drop_vols();
            
            %%%%%%%%%%%%
            %02_GradCorrect
            %For gradient non-linearity correction
            obj.Params.GradNonlinCorrect.in.movefiles = '../02_GradCorrect/';
            obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
            obj.Params.GradNonlinCorrect.in.gradfile = obj.gradfile;
            
            obj.Params.GradNonlinCorrect.in.fslroi = [ 0 1 ]; %To extraact the 1st b0
            obj.Params.GradNonlinCorrect.in.fn = obj.Params.DropVols.out.fn;

            obj.proc_gradient_nonlin_correct();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For FreeSurfer Segmentation (needed for b0 bbreg correction)
            [ tmpa, tmpb ] = system('whoami ');
            %[ tmpa, tmpb ] = system('echo $0');
            %[~,  tmpshell , ~] = fileparts(tmpb);
            obj.Params.FreeSurfer.shell = strtrim(tmpb); %strtrim(tmpshell);
            obj.Params.FreeSurfer.dir = obj.FS_location;
            obj.Params.FreeSurfer.init_location = obj.init_FS;
            
        
            %Retrieving a T1 scan:
            [sys_error, obj.Params.FreeSurfer.in.T1raw ] = system(['ls ' obj.session_location 'T1' filesep '*1mm.nii | head -1' ]);
            obj.Params.FreeSurfer.in.T1raw = strtrim(obj.Params.FreeSurfer.in.T1raw );
            if sys_error ~= 0 %No problem, we get the T1 the continue...
                fprintf(['\nError when finding the T1:'  obj.Params.FreeSurfer.in.T1raw  '\n'])
            end
            
            %Retrieving a T2 scan:
            [sys_error, obj.Params.FreeSurfer.in.T2_tempraw ] = system(['ls ' obj.session_location 'other' filesep '*T2SPACE* | head -1' ]);
            obj.Params.FreeSurfer.in.T2_tempraw = strtrim(obj.Params.FreeSurfer.in.T2_tempraw);
            if sys_error ~= 0 %No problem, we get the T1 the continue...
                fprintf(['\nNo T2 found:'  obj.Params.FreeSurfer.in.T2_tempraw  '\n'])
                obj.Params.FreeSurfer.in.T2_temprawexist=false;
            else
                [obj.Params.FreeSurfer.in.T1_dir, ~, ~ ] = fileparts(obj.Params.FreeSurfer.in.T1raw);
                obj.Params.FreeSurfer.in.T2exist=true;
                obj.Params.FreeSurfer.in.T2_dir = strrep(obj.Params.FreeSurfer.in.T1_dir,'T1','T2');
                system(['mkdir -p ' obj.Params.FreeSurfer.in.T2_dir ]);
                obj.Params.FreeSurfer.in.T2raw = [ obj.Params.FreeSurfer.in.T2_dir filesep 'T2.nii'];
                system(['cp '  strtrim(obj.Params.FreeSurfer.in.T2_tempraw) ' ' strtrim(obj.Params.FreeSurfer.in.T2raw)]);
            end
            
            obj.Params.FreeSurfer.out.aparcaseg = [ obj.Params.FreeSurfer.dir ...
                filesep obj.sessionname filesep 'mri' filesep 'aparc+aseg.mgz' ] ;
            
            obj.proc_getFreeSurfer();
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Extracting the value sfrom FreeSurfer:
            obj.getdata_FreeSurfer();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For b0 motion correciton (based on interspersed b0s)
            obj.fsdir=[obj.FS_location obj.sessionname ] ; 
            obj.Params.B0MoCo.FS = obj.fsdir;
            obj.Params.B0MoCo.in.movefiles = '../03_B0s_MoCo/';
            obj.Params.B0MoCo.in.prefix = 'moco_';
            obj.Params.B0MoCo.in.nDoF = '12' ; 
            obj.Params.B0MoCo.in.grad_rel = obj.Params.GradNonlinCorrect.out.warpfile;
            
            obj.Params.B0MoCo.in.fn = obj.Params.GradNonlinCorrect.out.fn;
            obj.Params.B0MoCo.in.bvals = obj.Params.DropVols.out.bvals;
            obj.Params.B0MoCo.in.bvecs = obj.Params.DropVols.out.bvecs;
            obj.Params.B0MoCo.in.sh_rotate_bvecs = obj.b0MoCo_rotate_bvecs_sh; 
            obj.proc_b0s_MoCo();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For BET2:
            obj.Params.Bet2.in.movefiles = ['..' filesep '04_Bet'];
            obj.Params.Bet2.in.fracthrsh = 0.4;
            obj.Params.Bet2.in.fn = obj.Params.B0MoCo.out.fn;
            
            obj.proc_bet2();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For EDDY:
            obj.Params.Eddy.in.movefiles = ['..' filesep '04_Eddy'];
            obj.Params.Eddy.in.fn=obj.Params.B0MoCo.out.fn;
            obj.Params.Eddy.in.bvals=obj.Params.B0MoCo.out.bvals';
            obj.Params.Eddy.in.bvecs=obj.Params.B0MoCo.out.bvecs;
            obj.Params.Eddy.in.mask = obj.Params.Bet2.out.mask;
            obj.Params.Eddy.in.index= ones(1,67); %for 35 volumes
            obj.Params.Eddy.in.acqp= [ 0 -1 0 0.08201 ]; %PE=A>>P (-1 at Y ) Echo spacing = 0.59 and EPI factor = 140 ==> 0.59(EPI)*0.001*(PE-1)139
            
            obj.proc_eddy();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %DERIVED MOVEMENT FROM EDDY:
            obj.Params.EddyMotion.in.movefiles = ['..' filesep '05_MotionFromEDDY'];
            obj.Params.EddyMotion.in.fn_eddy = obj.Params.Eddy.out.fn ;
            
            obj.proc_get_eddymotion(); obj.resave();
                        
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %To generate a mask after Eddy:
            %This step will 1) define a better mask if eddy affected the
            %movement of the head and 2) remove issues known to happen at
            %the edges of the brain when using the --wls option in dtifit!\
            %Reference: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;6eb4d787.1610
            obj.Params.MaskAfterEddy.in.movefiles = ['..' filesep '05_MaskAfterEddy'];
            obj.Params.MaskAfterEddy.in.fn = obj.Params.Eddy.in.fn; %Since we don't have a b0, we pass the full dwi and the method will take care of
            obj.Params.MaskAfterEddy.in.prefix = 'after_eddy';
            obj.Params.MaskAfterEddy.in.fracthrsh = '0.4';
            
            obj.proc_mask_after_eddy();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Until now, we have treated all images separately, so now we
            %will rigid-body coregister all *using flirst and spline
            %interpolation (proved to be more accurate):
            obj.Params.CoRegMultiple.in.fn = obj.Params.Eddy.in.fn;
            obj.Params.CoRegMultiple.in.b0 = obj.Params.MaskAfterEddy.in.b0 ;
            obj.Params.CoRegMultiple.in.bvals = obj.Params.Eddy.in.bvals;
            obj.Params.CoRegMultiple.in.bvecs = obj.Params.Eddy.out.bvecs;
            
            obj.Params.CoRegMultiple.in.movefiles = ['..' filesep '06_CoRegDWIs'];
            obj.Params.CoRegMultiple.in.ref_iteration = 2; % All images will be registered to this iteration (in ADRC, 7p5_set1, index 1 is for 2p7_set4!)
            
            
            obj.proc_coreg_multiple();
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For DWI combined DTIFIT:
            %We will use the --wls option as it seems to improve the fit of
            %the diffusion tensor model and negativity values due to noise
            %REF: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;735f6320.1409
            %ALSO WE WILL ONLY USE THE 2.5k because w/ 7.5k the values
            %diverge for linearity (2.5k is close to linearity than
            %non-linear signal drop at 7.5k)
            obj.Params.Dtifit.in.movefiles = [ '..' filesep '06_2k5DTIFIT' ];
            obj.Params.Dtifit.in.fn ={obj.Params.CoRegMultiple.out.fn{1}}; %{obj.Params.CoRegMultiple.out.combined_fn};
            obj.Params.Dtifit.in.prefix = 'DTIFIT_FSLv509' ; %Double check this so you prefix the version of FSL!
            obj.Params.Dtifit.in.bvecs = {obj.Params.CoRegMultiple.out.bvecs{1}};
            obj.Params.Dtifit.in.bvals = {obj.Params.CoRegMultiple.out.bvals{1}};
            obj.Params.Dtifit.in.mask = {obj.Params.CoRegMultiple.out.combined_mask};
            
            obj.proc_dtifit();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For GQI:
            obj.Params.GQI.in.movefiles = [ '..' filesep '07_Combined_Recon_gqi' ];
            obj.Params.GQI.in.fn = {obj.Params.CoRegMultiple.out.combined_fn};
            obj.Params.GQI.in.bvecs = {obj.Params.CoRegMultiple.out.combined_bvecs};
            obj.Params.GQI.in.bvals = {obj.Params.CoRegMultiple.out.combined_bvals};
            obj.Params.GQI.in.mask ={obj.Params.CoRegMultiple.out.combined_mask};
            obj.Params.GQI.in.prefix = 'GQI_DSISv053117' ; %Double check this so you prefix the version of DSISTUDIO!
            obj.Params.GQI.out.export = 'gfa,nqa0,nqa1';
            
            obj.proc_gqi();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %FS2dwi:
            obj.Params.FS2dwi.in.movefiles = ['..' filesep '07_Combined_FS2dwi' ];
            obj.Params.FS2dwi.in.b0 =  {obj.Params.CoRegMultiple.out.combined_bet} ; % Removed due to nonskulls stripping --> {obj.Params.CoRegMultiple.out.combined_b0} ; 
            obj.Params.FS2dwi.in.aparcaseg = obj.Params.FreeSurfer.out.aparcaseg ; 
            
            obj.Params.FS2dwi.in.tmpfile_aparcaseg = [ obj.dependencies_dir filesep 'FS_DEPS' filesep  'FS_aparc.txt' ] ; 
            obj.Params.FS2dwi.in.tmpfile_aparcaseg2009 = [ obj.dependencies_dir filesep 'FS_DEPS' filesep   'FS_aparc2009.txt' ] ; 
            obj.Params.FS2dwi.in.tmpfile_hippo_bil = [ obj.dependencies_dir filesep 'FS_DEPS' filesep   'FS_hippolabels_bil.txt' ] ;
            
            
            
            obj.Params.FS2dwi.in.aparcaseg2009 = ...
                strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','aparc.a2009s+aseg')); 
            
            %A possible error is the naming convention when only a T1 was
            %used!!
            obj.Params.FS2dwi.in.hippofield_left = ...
                strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','lh.hippoSfLabels-T1-T2.v10')); 
            obj.Params.FS2dwi.in.hippofield_right = ...
                strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','rh.hippoSfLabels-T1-T2.v10')); 
          
            obj.proc_FS2dwi();
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %FROIS2dwi:
            obj.Params.FROIS2dwi.in.movefiles = [ '..' filesep '07_FROIS2dwi' ];
            obj.Params.FROIS2dwi.in.fn = obj.Params.CoRegMultiple.out.combined_bet;
            obj.Params.FROIS2dwi.in.prefix = 'MNIT1_2_dwi' ; %Double check this so you prefix the version of DSISTUDIO!
            obj.Params.FROIS2dwi.in.FROIS_dir = obj.FROIS_dir;
            obj.Params.FROIS2dwi.in.MNI_T1 = [ obj.FROIS_dir 'template' filesep 'MNI152_T1_1mm_brain.nii.gz' ] ;
          
           %On the works.... 
           %obj.proc_FROIS2dwi();
        end
  
        function obj = CommonPostProc(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Creating Fornix TRKLAND
            %~~~~~~~ TRACULA (and implicit functionality of bedpostx):
            for tohide=1:1
            obj.Params.Tracula.in.movefiles = ['..' filesep 'post_TRACULA' ];
            obj.Params.Tracula.in.fn = obj.Params.CoRegMultiple.out.combined_fn ;
            obj.Params.Tracula.in.dcmrirc = [obj.dependencies_dir filesep 'TRACULA_DEPS' filesep 'dcmrirc.template' ];
            obj.Params.Tracula.in.FSDIR = obj.Params.FreeSurfer.dir;
            obj.Params.Tracula.in.bvec = obj.Params.CoRegMultiple.out.combined_bvecs ;
            obj.Params.Tracula.in.bval = obj.Params.CoRegMultiple.out.combined_bvals ;
            obj.Params.Tracula.in.nb0 = 28;
            obj.Params.Tracula.in.prefix = 'adrc';
            
            obj.proc_tracula(); obj.resave();
            end
           
            
            %~~~~~~~ White matter Lesions 2 DWI space:
            for tohide=1:1
                obj.Params.WMLs2DWI.in.movefiles = ['..' filesep 'post_WML2DWI' ];
                obj.Params.WMLs2DWI.in.b0 = obj.Params.CoRegMultiple.out.combined_b0;
                
                obj.Params.WMLs2DWI.in.dir = '/eris/bang/ADRC/PROJECTS/WMLs_LST_LGA/FLAIRS/' ;
                obj.Params.WMLs2DWI.in.FLAIR = [obj.Params.WMLs2DWI.in.dir 'm' obj.sessionname '_FLAIR.nii' ];
                obj.Params.WMLs2DWI.in.WMLprobmap = [obj.Params.WMLs2DWI.in.dir 'ples_lpa_m' ...
                    obj.sessionname '_FLAIR.nii' ];
                obj.proc_WMLs2DWI(); obj.resave();
            end
        
            %~~~~~~ TRKLAND PROCESSING:
            obj.Trkland.root = [ obj.root  'post_TRKLAND' filesep ];
            % TRKLAND_FX:
            for tohide=1:1
                obj.Trkland.fx.in.movefiles = ['..' filesep 'post_TRKLAND' ];
                %b0 params:
                %obj.Trkland.fx.in.b0 = obj.Params.CoRegMultiple.out.combined_bet; %this will not be used anymore for coreg as FA proves to do a more optimal job (below)
                obj.Trkland.fx.in.b0 = obj.Params.CoRegMultiple.out.combined_b0;
                obj.Trkland.fx.in.FA = obj.Params.Dtifit.out.FA{1};
                %Temobj.Params.CoRegMultiple.out.combined_b0plate parameters:
              
                %Based on orientation, we will chooose a specific template
                %(usually RAS) 
                [~, Ori ] = system(['mri_info ' obj.Trkland.fx.in.b0 ' | grep Orientation | awk ''{print $3}'''] );
                obj.Trkland.fx.tmp.ori = strtrim(Ori);
                clear Ori
                if strcmp(obj.Trkland.fx.tmp.ori,'LPS')
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
                obj.trkland_fx(); obj.resave();
            end
            
            % [TO MODIFY ROIs/ROAs] TRKLAND_HIPPOCING:
            for tohide=1:1
                obj.Trkland.hippocing.in.hippo_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Left-Hippocampus.nii.gz');
                obj.Trkland.hippocing.in.hippo_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Right-Hippocampus.nii.gz');
                obj.Trkland.hippocing.in.postcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-posteriorcingulate.nii.gz');
                obj.Trkland.hippocing.in.postcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-posteriorcingulate.nii.gz');
                
                %Interpolation n (for cingulum):
                obj.Trkland.hippocing.in.n_interp=33;
                obj.trkland_hippocing(); obj.resave();
            end
            
            %  [TO MODIFY ROIs/ROAs] TRKLAND_CINGULUM:
            for tohide=1:1
                obj.Trkland.cingulum.in.rostantcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-rostralanteriorcingulate.nii.gz');
                obj.Trkland.cingulum.in.rostantcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-rostralanteriorcingulate.nii.gz');
                obj.Trkland.cingulum.in.postcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-posteriorcingulate.nii.gz');
                obj.Trkland.cingulum.in.postcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-posteriorcingulate.nii.gz');
                
                % BASE ON THESE VALUES (FROM EARLIER ADRC_PROCESSING W/O INTERP:
                % ninter_fx = 40;
                % ninter_cingulum = 32;
                % ninter_hippocing = 33;
                obj.Trkland.cingulum.in.n_interp = 32;
                trkland_cingulum(obj); obj.resave();
            end
            
            %% [ DEPRECATED METHODS ]
            %~~~~~~ TRKLAND_DEPENDENT:
            % [ DEPRECATED ] TRKLAND_ATR:
             for tohide=1:1
%                 CANNOT USE TRKLAND_ATR APPROACH FOR THE ATR TRACTS AS
%                 A MOST ROBUST STREAMLINE CANNOT BE ISOLATED
%                 obj.Trkland.atr.in.rostantcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-rostralanteriorcingulate.nii.gz');
%                 obj.Trkland.atr.in.rostantcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-rostralanteriorcingulate.nii.gz');
%                 obj.Trkland.atr.in.thalamus_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Left-Thalamus-Proper.nii.gz');
%                 obj.Trkland.atr.in.thalamus_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Right-Thalamus-Proper.nii.gz');
%                 
%               %Under DEVELOPMENT:
%               %I haven't found an optimal trimming technique...yet
%               %trkland_atr(obj)
             end
             
            %~~~~~~ TRACULA DEPENDENT:
            % [ DEPRECATED ] TRACX THAL_2_CORTEX11:
            for tohide=1:1
%                 obj.Params.tracx_thal2ctx11.in.bedp_dir = fileparts(obj.Params.Tracula.out.bedp_check);
%                 obj.Params.tracx_thal2ctx11.in.FSaparc_dir = [ fileparts(obj.Params.FS2dwi.out.fn_aparc2009)  filesep 'aparc_aseg' filesep];
%                 obj.Params.tracx_thal2ctx11.in.movefiles = ['..' filesep '..' filesep '..' filesep 'post_tracx' filesep 'thal2ctx11' ];
%                 obj.Params.tracx_thal2ctx11.in.prep_segs_list = [ obj.dependencies_dir filesep 'TRACULA_DEPS' filesep  'THALX_CTX11.txt' ];
%               
               %obj.proc_tracx2thal11();
            end
            % [DEPRECATED ] TRACX THAL_2_PAPEZ (2 frontals, 1 cingualte, 3 temporals):
            for tohide=1:1
%                 obj.Params.tracx_thal2papez.in.bedp_dir = fileparts(obj.Params.Tracula.out.bedp_check);
%                 obj.Params.tracx_thal2papez.in.FSaparc_dir = [ fileparts(obj.Params.FS2dwi.out.fn_aparc2009)  filesep 'aparc_aseg' filesep];
%                 obj.Params.tracx_thal2papez.in.movefiles = ['..' filesep '..' filesep '..' filesep 'post_tracx' filesep 'thal2papez' ];
%                 obj.Params.tracx_thal2papez.in.prep_segs_list = [ obj.dependencies_dir  filesep 'TRACULA_DEPS' filesep 'ontheworks' filesep  'THALX_PAPEZ.txt' ];
%                 
%                 obj.proc_tracx2papez();
            end
            % [ DEPRECATED ] TRACX THAL_2_DMN:
            for tohide=1:1
%                 obj.Params.tracx_thal2dmn.in.bedp_dir = fileparts(obj.Params.Tracula.out.bedp_check);
%                 obj.Params.tracx_thal2dmn.in.DMN_dir = obj.FROIS_dir ; 
%                 obj.Params.tracx_thal2dmn.in.movefiles = ['..' filesep '..' filesep '..' filesep 'post_tracx' filesep 'thal2ctx10' ];
%                 obj.Params.tracx_thal2dmn.in.prep_segs_list = [ obj.dependencies_dir  filesep 'TRACULA_DEPS' filesep  'THALX_CTX10.txt' ];  
%                 
%                 obj.proc_tracx2DMN();
            end
            
            %~~~~~~~ AFQ (https://github.com/yeatmanlab/AFQ/wiki) ADAPTABILITY DEPRECATED 
            %(NO MANUAL AC-PC WHICH IS NECESSARY FOR AFQ TO WORK)
            % [DEPRECATED ] AFQ:
            for tohide=1:1
                % obj.Params.AFQ.in.movefiles = ['..' filesep 'post_AFQ' ];
                % obj.Params.AFQ.in.dwi = obj.Params.CoRegMultiple.out.combined_fn ;
                % obj.Params.AFQ.in.T1 = obj.T1 ;
                %
                % obj.proc_AFQ();
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
      
    end
    methods ( Access = protected )
        function obj = getDCM2nii(obj,torun)
            %For proc_DCM2NII:
            obj.Params.DCM2NII.specific_vols=68;
            obj.Params.DCM2NII.scanlog = [ obj.session_location filesep 'LogFiles' ...
                filesep 'scan.log' ] ;
            if ~exist(obj.Params.DCM2NII.scanlog,'file')
                error(['No scanlog file found in: ' obj.Params.DCM2NII.scanlog ' . Exiting...']);
            end
            obj.Params.DCM2NII.seq_names={ 'ep2d_diff_7p5k_set1E60' 'ep2d_diff_7p5k_set2E60' ...
                'ep2d_diff_7p5k_set3E60' 'ep2d_diff_2p5k_set4E60' };
            
           
            for ii=1:4 % 4 sets of DWIs in this project!
                  obj.Params.DCM2NII.in(ii).fsl2std_param = '-1 0 0 250.199 \n0 1 0 250.199 \n0 0 -1 0 \n0 0 0 1';
           
                obj.Params.DCM2NII.in(ii).prefix = obj.Params.DCM2NII.seq_names(ii);
                [ ~ , obj.Params.DCM2NII.in(ii).nvols ] = system([ 'cat ' ...
                    obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names{ii} ...
                    ' | tail -1 | awk ''{ print $7 }'' ' ]);
                obj.Params.DCM2NII.in(ii).nvols=str2num(obj.Params.DCM2NII.in(ii).nvols);
                [ ~ , obj.Params.DCM2NII.in(ii).first_dcmfiles ]= system([ 'cat ' ...
                    obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names{ii} ...
                    ' | tail -1 | awk ''{ print $8 }'' ' ]);
                
                obj.Params.DCM2NII.out(ii).location = [ obj.root 'Orig' filesep ];
                obj.Params.DCM2NII.out(ii).fn = [  obj.Params.DCM2NII.out(ii).location  cell2char(obj.Params.DCM2NII.seq_names(ii)) '.nii.gz' ];
            end
            if (torun) ; obj.proc_dcm2nii ; end 
        end
    end
end