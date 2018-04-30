classdef dwi_HABSIIC < dwiMRI_Session
    %%  This class is a subclass of its parent class dwi_MRI_Session.m
    %%  (where it will inherent other methods).
    %%  Created by: RodrigoPereaCamargo ~~> rperea14@live.com
    %%  Inheritance credits to AaronSchultz
    %%

    properties
       %Properties initialized in parent class dwiMRI_Session() and 
        %And assgiend values in method: obj.setMyParams() which is called
        %in the class constructor dwi_ADRC(). 
    end
    
    methods
        function obj = dwi_HABSIIC(sessionname,opt)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %INITIALIZE VARIABLES NOW:
            obj.setMyParams();
            
            %Initialize root variables:
            obj.projectID='HABSIIC';
            obj.sessionname = sessionname;
            obj.root = [obj.session_dir sessionname '/DWIs/'];
            obj.dcm_location = [ obj.dcm_location sessionname filesep ];
            obj.session_location= [ obj.session_dir sessionname filesep ] ;
           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Check to see if a *.mat file exists.
            if exist([obj.session_location 'DWIs' filesep  sessionname '.mat'],'file')>0
                load([obj.session_location 'DWIs' filesep  sessionname '.mat']);
                oldroot = obj.root;
                obj.wasLoaded = true;
            else
                %  obj.setMyParams;
                donothing=1;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Check if you are in the right project (if not, probably obj.session_location doesn't exist:
            if ~exist([obj.session_location ],'dir')
                error(['Sesion directory: ' obj.session_location sessionname ' doesn''t exist. Are you sure are in the right Project ID?']);
            end

            
            %%%%%%%%%%%%%%%%%%%%%%%%
            %Utilizing 1st Argument:
            if nargin>1
                obj.root = opt;
            end
             %If the folder /DWIs/ does not exist, then create it!
            if exist(obj.root,'dir')==0
                obj.make_root();
            end
            obj.objectHome = obj.root ;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %CHECK CHANGES MADE FROM /eris to /cluster
            %CODE CAN BE RECYCLE TO CHANGE VARIABLE STRINGS WITHIN THE
            %OBJECT
%             if ~strcmp(strtrim(obj.FS_location),'/autofs/cluster/hab/HAB/HABS_IIC/Fresurfer/')
%                 display('Changing eris to cluster folder...');
%                 %Then apply the mod to cluster...
%                 newobj = replaceObjText_v2(obj,{'eris'},{'cluster'});
%                 obj=newobj;
%                 obj.objectHome = obj.root ;
%                 display('CHANGING: ''eris'' was changed to ''cluster'' ');
%                 obj.resave();
%             end
           
            %Check if *.nii.gz files exist, if not get them from DCM2nii:
%             obj.rawfiles = dir_wfp([obj.root 'Orig/*.nii.gz' ] );
%             if isempty(obj.rawfiles) || numel(obj.rawfiles) ~= 4 % 4 DWIs sequence acquired here
%                 RunFlag = true ;
%                 obj.getDCM2nii(RunFlag);
%             end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %Start the CommonPreProc:
%             obj.CommonPreProc();
%             %Start the CommonPostProc (commenting it? Maybe, maybe not?):
%             obj.CommonPostProc();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             if isdeployed()
%                 display('Diffusion Pre- and Post- processing completed')
%                 fprintf('\n\n *If running the compiled version and for completion, \n make sure you run (from MATLAB terminal):')
%                 fprintf('\n * 1. Load the object ~~> load(<Sessions>/<SUBJID>/DWIs/<SUBJID>.mat) \n')
%                 fprintf('\n * 2. Run             ~~> obj.proc_T1toDWI() \n')
%                 fprintf('\n * 3. Run             ~~> obj.prep_spmT1_proc() \n')
%                 fprintf('\n * 4. Run             ~~> obj.trkland_fx() \n')
%             end
            
        end
        
        function obj= setMyParams(obj)
            %Remove previous variable no needed anymore:
            obj.T1 = 'Please refer to --> obj.Params.T1toDWI.in.T1';
            
            %Check protected method initProperties()
            %root directoy where raw data lives:
            obj.object_dir= '/cluster/brutha/MATLAB_Scripts/Objects/Diffusion/'; %To add to path if needed
            obj.session_dir='/cluster/hab/HAB/HABS_IIC/Sessions/';
            obj.dcm_location = '/cluster/hab/HAB/HABS_IIC/DICOM_ARCHIVE/';
            obj.gradfile='/autofs/cluster/hab/HAB/HABS_IIC/Scripts/coeff_AS82.grad';
            obj.dependencies_dir='/cluster/bang/ADRC/Scripts/DEPENDENCIES/'; % Utility Scripts
            
            %FreeSurfer Dependencies
            obj.FS_location='/cluster/hab/HAB/HABS_IIC/Fresurfer/';
            obj.init_FS = '/usr/local/freesurfer/stable6/';
            obj.FSL_dir = '/usr/pubsw/packages/fsl/5.0.9/';
            
            %Trkland dependencies:
            obj.fx_template_dir='/autofs/cluster/bang/ADRC/TEMPLATES/FX_1.8mm_orig/';
            
            %Dependency properties:
%             obj.sh_gradfile=['/cluster/bang/ADRC/Scripts/DEPENDENCIES/GradNonLin_Correc/run_mris_gradient_nonlin__unwarp_volume__batchmode_ADRC_v3.sh ' ...
%                 '/usr/pubsw/common/matlab/8.5'];
            obj.b0MoCo_rotate_bvecs_sh='/cluster/bang/ADRC/Scripts/DEPENDENCIES/PREPROC_DEPS/rotate_bvecs.sh'; %For rotating the bvecs after proc_b0MoCo
            
            obj.init_rotate_bvecs_sh='/cluster/bang/ADRC/Scripts/DEPENDENCIES/PREPROC_DEPS/mod_fdt_rotate_bvecs.sh'; %THIS IS ONLY USED IN THE proc_dcm2nii() METHOD!
            obj.col2rows_sh='/cluster/bang/ADRC/Scripts/DEPENDENCIES/PREPROC_DEPS/drigo_col2rows.sh';
            obj.redo_history = false; %Allows to redo the history of all processes withouth running any obj.BashCode. IT SHOULD ALWAYS BE FALSE UNLESS OTHERWISE!
        end
        
        function resave(obj)
            if ~(obj.redo_history)
                save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
            end
        end
        
        %COMMON METHODS THAT CALL post_XXX METHODS IN SUPERCLASS dwiMRI_Session.m:
        function obj = CommonPreProc(obj)
            %keyboard; 
            if isempty(obj.rawfiles)
                obj.rawfiles = dir_wfp([obj.root 'Orig' filesep '*.nii*' ] );
                if isempty(obj.rawfiles)
                    error('Cannot find any *.nii* files in obj.root/Orig/*.nii*. Please check')
                end
            end           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %01_Merging all the AP and PAs DWIs
            fprintf('MERGING DWIs...\n');
            obj.Params.mergeDWIs.root =[obj.root '01_MergedDWIs/'];
            if exist(obj.Params.mergeDWIs.root,'dir') == 0 
                mkdir(obj.Params.mergeDWIs.root);
            end
            %AP:
            for tohide=1:1
                % INPUT:
                obj.Params.mergeDWIs.AP.in.dwi = dir_wfp([obj.root 'Orig/dMRI*AP*.nii*']) ;
                for ii=1:numel(obj.Params.mergeDWIs.AP.in.dwi)
                    obj.Params.mergeDWIs.AP.in.bval{ii} = strrep(obj.Params.mergeDWIs.AP.in.dwi{ii},'.nii','.bval');
                    obj.Params.mergeDWIs.AP.in.bvec{ii} = strrep(obj.Params.mergeDWIs.AP.in.dwi{ii},'.nii','.bvec');
                end
                % OUTPUT:
                obj.Params.mergeDWIs.AP.out.dwi  = [ obj.Params.mergeDWIs.root 'merged_dMRI_AP.nii' ];
                obj.Params.mergeDWIs.AP.out.bval = [ obj.Params.mergeDWIs.root 'merged_dMRI_AP.bvals' ];
                obj.Params.mergeDWIs.AP.out.bvec = [ obj.Params.mergeDWIs.root 'merged_dMRI_AP.voxel_space.bvecs' ];
            end
            %Run the common method to merge:
            obj.proc_mergeDWIs(obj.Params.mergeDWIs.AP.in,obj.Params.mergeDWIs.AP.out);
            
            %PA:
            for tohide=1:1
                % INPUT:
                obj.Params.mergeDWIs.PA.in.dwi = dir_wfp([obj.root 'Orig/dMRI*PA*.nii*']) ;
                for ii=1:numel(obj.Params.mergeDWIs.PA.in.dwi)
                    obj.Params.mergeDWIs.PA.in.bval{ii} = strrep(obj.Params.mergeDWIs.PA.in.dwi{ii},'.nii','.bval');
                    obj.Params.mergeDWIs.PA.in.bvec{ii} = strrep(obj.Params.mergeDWIs.PA.in.dwi{ii},'.nii','.bvec');
                end
                % OUTPUT:
                obj.Params.mergeDWIs.PA.out.dwi  = [ obj.Params.mergeDWIs.root 'merged_dMRI_PA.nii' ];
                obj.Params.mergeDWIs.PA.out.bval = [ obj.Params.mergeDWIs.root 'merged_dMRI_PA.bvals' ];
                obj.Params.mergeDWIs.PA.out.bvec = [ obj.Params.mergeDWIs.root 'merged_dMRI_PA.voxel_space.bvecs' ];
            end
            %Run the common method to merge:
            obj.proc_mergeDWIs(obj.Params.mergeDWIs.PA.in,obj.Params.mergeDWIs.PA.out);
            fprintf('DONE MERGING DWIs...\n');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %02_GNC
            %% check that the AP/PA files are chosen correctly.
            if isempty(dir_wfp([obj.root '02_GNC/gnc_merged_dMRI*AP.nii']))
                obj.Params.GradNonlinCorrect.in.movefiles = '';
                obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
                obj.Params.GradNonlinCorrect.in.gradfile = '/autofs/cluster/hab/HAB/HABS_IIC/Scripts/coeff_AS82.grad';
                fn = dir_wfp([obj.root '01_MergedDWIs/merged_dMRI*AP*.nii*']);
                obj.Params.GradNonlinCorrect.in.fn = fn;
                if exist([obj.root '02_GNC/AP_dwi_grads'],'dir') == 0
                    mkdir([obj.root '02_GNC/AP_dwi_grads']);
                end
                [m, h] = read_image([fn{1} ',1'],0);
                h = h(1);
                h.fname = [obj.root '02_GNC/AP_dwi_grads/firstb0.nii'];
                spm_write_vol(h,m);
                
                obj.Params.GradNonlinCorrect.in.target = h.fname;

                obj.Params.GradNonlinCorrect.out.warpfile = [];
                obj.Params.GradNonlinCorrect.out.meannii = [];
                obj.Params.GradNonlinCorrect.out.fn = [];
                obj.proc_gradient_nonlin_correct();
                
                system(['mv ' obj.root '01_MergedDWIs/gnc* ' obj.root '01_GNC/']);
            end
            %%
             if isempty(dir_wfp([obj.root '02_GNC/gnc_merged_dMRI*PA.nii']))
            
                obj.Params.GradNonlinCorrect.in.movefiles = '';
                obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
                obj.Params.GradNonlinCorrect.in.gradfile = '/autofs/cluster/hab/HAB/HABS_IIC/Scripts/coeff_AS82.grad';
                fn = dir_wfp([obj.root '01_MergedDWIs/merged_dMRI*PA*.nii*']);
                obj.Params.GradNonlinCorrect.in.fn = fn;
                
                mkdir([obj.root '02_GNC/PA_dwi_grads']);
                
                [m, h] = read_image([fn{1} ',1'],0);
                h = h(1);
                h.fname = [obj.root '02_GNC/PA_dwi_grads/firstb0.nii'];
                spm_write_vol(h,m);
                
                obj.Params.GradNonlinCorrect.in.target = h.fname;

                obj.Params.GradNonlinCorrect.out.warpfile = [];
                obj.Params.GradNonlinCorrect.out.meannii = [];
                obj.Params.GradNonlinCorrect.out.fn = [];
                obj.proc_gradient_nonlin_correct();
                
                system(['mv ' obj.root '01_MergedDWIs/gnc* ' obj.root '02_GNC/']);
            end
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %03_B0 motion correciton (based on interspersed b0s)
            obj.Params.B0MoCo.in.movefiles = '../03_B0s_MoCo/';
            obj.Params.B0MoCo.in.fn = dir_wfp([obj.root '02_GNC/gnc_merged_dMRI*.nii']);
            %Assigning bvecs, bvals:
            bvals = []; bvecs = [];
            for ii = 1:numel(obj.Params.B0MoCo.in.fn)
                [s1, s2, s3] = fileparts(obj.Params.B0MoCo.in.fn{ii});
                obj.Params.B0MoCo.in.bvals{ii,1} = ([obj.root '01_MergedDWIs/' regexprep(s2,'gnc_','') '.bvals']);
                obj.Params.B0MoCo.in.bvecs{ii,1} = ([obj.root '01_MergedDWIs/' regexprep(s2,'gnc_','') '.voxel_space.bvecs']);
            end
            obj.Params.B0MoCo.in.prefix = 'moco_';
            
            obj.Params.B0MoCo.in.refB0 = 1 ; %First b0! 
            obj.Params.B0MoCo.in.sh_rotate_bvecs = obj.b0MoCo_rotate_bvecs_sh;
            
            obj.proc_b0s_MoCo(1); %First argument '1' denotes the refB0 is the first volume (in index FSL 0000.nii.gz)
            %%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %04_B0mean:
            obj.Params.B0mean.in.movefiles = ['..' filesep '04_B0mean'];
            obj.Params.B0mean.in.fn= obj.Params.B0MoCo.out.fn;
            
            obj.Params.B0mean.in.b0_nvols=[]; %states number of volumes at the beginning. Is empty, select b0 volumes automatically
            %instead use the location of bvals:
            obj.Params.B0mean.in.fn_bvals=obj.Params.B0MoCo.out.bvals;
            %and a b0indx:
            obj.Params.B0mean.in.b0idx =5;
            obj.proc_meanb0();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %05_TOPUP:
            obj.Params.Topup.in.movefiles = ['..' filesep '05_Topup'];
            obj.Params.Topup.in.fn_AP= dir_wfp([obj.root '03_B0s_MoCo/moco_gnc_merged_dMRI*AP.nii*']);
            obj.Params.Topup.in.fn_PA= dir_wfp([obj.root '03_B0s_MoCo/moco_gnc_merged_dMRI*PA.nii*']);
            obj.Params.Topup.in.b0_AP= dir_wfp([obj.root '04_B0mean/b0mean*AP.nii*']);
            obj.Params.Topup.in.b0_PA= dir_wfp([obj.root '04_B0mean/b0mean*PA.nii*']);
            % Echo spacing:  0.69ms
            % EPI factor: 140
            % So acqp is 0.001*0.69*139 = 0.096;
            obj.Params.Topup.in.acqp_params = [ 0 1 0 0.096 ; 0 -1 0 0.096] ;
            
            obj.proc_topup();
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For EDDY:
            obj.Params.Eddy.in.movefiles = ['..' filesep '06_Eddy'];
            
            %Merginfiles before executing Eddy:
            in_eddy.dwi = obj.Params.B0MoCo.out.fn';
            in_eddy.bval = obj.Params.B0MoCo.out.bvals';
            in_eddy.bvec = obj.Params.B0MoCo.out.bvecs';
            out_eddy.dwi = [ obj.root '06_Eddy/merged_AP_PA_dMRI.nii.gz' ];
            out_eddy.bval = [ obj.root '06_Eddy/merged_AP_PA_dMRI.bvals' ];
            out_eddy.bvec = [ obj.root '06_Eddy/merged_AP_PA_dMRI.voxel_space.bvecs' ];
            
            obj.proc_mergeDWIs(in_eddy, out_eddy);
            %...finished merging AP_PAs
            
            
            obj.Params.Eddy.in.fn{1}=out_eddy.dwi;
            obj.Params.Eddy.in.bvals{1}=out_eddy.bval;
            obj.Params.Eddy.in.bvecs{1}=out_eddy.bvec;
            obj.Params.Eddy.in.mask{1} =obj.Params.Topup.out.bet_mask;
            obj.Params.Eddy.in.index= [ ones(1,199) 2*ones(1,199) ]; %for 35 volumes
            obj.Params.Eddy.in.topup = obj.Params.Topup.out.cnf_prefix ; 
            
            %RePLACED
            %obj.Params.Eddy.in.acqp= [ 0 -1 0 0.08201 ]; %PE=A>>P (-1 at Y ) Echo spacing = 0.59 and EPI factor = 140 ==> 0.59(EPI)*0.001*(PE-1)139
            %WITH (only in this porject that uses TOPUP, the others will
            %create it automatically give the .in.acqp...
            obj.Params.Eddy.out.fn_acqp{1} =[ obj.root '06_Eddy/acqp.txt' ];
            if exist(obj.Params.Eddy.out.fn_acqp{1} ,'file') == 0 
                %AP direction:
                for ii=1:199
                    system(['printf " 0 1 0 0.096\n" >> ' obj.Params.Eddy.out.fn_acqp{1}  ]);
                end
                %PA direction:
                for ii=1:199
                    system(['printf " 0 -1 0 0.096\n" >> ' obj.Params.Eddy.out.fn_acqp{1}]);
                end
            end
            obj.proc_eddy();
            keyboard
            %%
            %%%%%%%%%%%%
%             %02_GradCorrect
%             % For gradient non-linearity correction
%             obj.Params.GradNonlinCorrect.in.movefiles = '../02_GradCorrect/';
%             obj.Params.GradNonlinCorrect.in.fn = obj.Params.DropVols.out.fn;
%             obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
%             obj.Params.GradNonlinCorrect.in.gradfile = obj.gradfile;
%             obj.Params.GradNonlinCorrect.in.fslroi = [ 0 1 ]; %To extraact the 1st b0 and apply gnc only to it
%             
%             obj.proc_gradient_nonlin_correct();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %??_Proc_FreeSurfer
%             obj.Params.FreeSurfer.init_location = obj.init_FS;
%             obj.Params.FreeSurfer.dir = obj.FS_location;
%             obj.Params.FreeSurfer.out.aparcaseg = [ obj.Params.FreeSurfer.dir ...
%                 obj.sessionname filesep 'mri' filesep 'aparc+aseg.mgz' ];
%             %T1_01_For FreeSurfer Segmentation (needed for segmenting ROIs to diffusion space)
%             %Retrieving a T1 scan:
%             [~ , obj.Params.FreeSurfer.in.T1raw ] = system(['ls ' obj.session_location 'T1' filesep 'rawT1_HCPF_1mm_001.nii | head -1' ]);
%             obj.Params.FreeSurfer.in.T1raw = strtrim(obj.Params.FreeSurfer.in.T1raw );
%             if exist(obj.Params.FreeSurfer.in.T1raw,'file') == 0 %No problem, we get the T1 the continue...
%                 error(['\nT1 not found:'  obj.Params.FreeSurfer.in.T1raw  '\n'])
%             end
%             %Retrieving a T2 scan:
%             [~ , obj.Params.FreeSurfer.in.T2raw ] = system(['ls ' obj.session_location 'T2' filesep 'rawT2_HCPF_1mm* | head -1' ]);
%             obj.Params.FreeSurfer.in.T2raw = strtrim(obj.Params.FreeSurfer.in.T2raw);
%             if exist(obj.Params.FreeSurfer.in.T2raw,'file') == 0 %No problem, we get the T1 the continue...
%                 fprintf(['\nT2 not found:'  obj.Params.FreeSurfer.in.T2raw  '\n'])
%                 obj.Params.FreeSurfer.in.T2exist=false;
%             else
%                 obj.Params.FreeSurfer.in.T2exist=true;
%             end
%             
%             %Some optional parameters specific to what SHELL envinroment is
%             %being used. If 'rdp20' then syntanx should be in bash.
%             %Otherwise, the syntax should be in csh/tsch
%             [ ~, tmpb ] = system('whoami ');
%             obj.Params.FreeSurfer.shell = strtrim(tmpb);
%             
%             obj.proc_getFreeSurfer();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %T1_02_Getting the necessary values from FreeSurfer's output:
            obj.getdata_FreeSurfer();
                        
            
            
            
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
            obj.Params.EddyMotion.in.fn_eddy = obj.Params.Eddy.out.fn;
            
            obj.proc_get_eddymotion();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %To generate a mask after Eddy:
            %This step will 1) define a better mask if eddy affected the
            %movement of the head and 2) remove issues known to happen at
            %the edges of the brain when using the --wls option in dtifit!\
            %Reference: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;6eb4d787.1610
            obj.Params.MaskAfterEddy.in.movefiles = ['..' filesep '05_MaskAfterEddy'];
            obj.Params.MaskAfterEddy.in.fn = obj.Params.Eddy.out.fn; %Since we don't have a b0, we pass the full dwi and the method will take care of the rest
            obj.Params.MaskAfterEddy.in.prefix = 'after_eddy';
            obj.Params.MaskAfterEddy.in.fracthrsh = '0.4';
            
            obj.proc_mask_after_eddy();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Until now, we have treated all images separately, so now we
            %will rigid-body coregister all *using flirst and spline
            %interpolation (proved to be more accurate):
            obj.Params.CoRegMultiple.in.movefiles = ['..' filesep '06_CoRegDWIs'];
            obj.Params.CoRegMultiple.in.fn = obj.Params.Eddy.in.fn;
            obj.Params.CoRegMultiple.in.b0 = obj.Params.MaskAfterEddy.in.b0;
            obj.Params.CoRegMultiple.in.bvals = obj.Params.Eddy.in.bvals;
            obj.Params.CoRegMultiple.in.bvecs = obj.Params.Eddy.out.bvecs;
            obj.col2rows_sh='/cluster/bang/ADRC/Scripts/DEPENDENCIES/PREPROC_DEPS/drigo_col2rows.sh';
            obj.b0MoCo_rotate_bvecs_sh ='/cluster/bang/ADRC/Scripts/DEPENDENCIES/PREPROC_DEPS/rotate_bvecs.sh';
            obj.col2rows_sh='/cluster/bang/ADRC/Scripts/DEPENDENCIES/PREPROC_DEPS/drigo_col2rows.sh';
            obj.Params.CoRegMultiple.in.ref_iteration = 2; % All images will be registered to this iteration (in ADRC, 7p5_set1, index 1 is for 2p7_set4!)
            
            
            obj.proc_coreg_multiple();
                        
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Coregistering the T1 to B0:
            obj.Params.T1toDWI.in.movefiles = './07_T1toDWI/';
            obj.Params.T1toDWI.in.b0 = obj.Params.CoRegMultiple.out.combined_b0;
            obj.Params.T1toDWI.in.T1 = obj.Params.FreeSurfer.in.T1;
           
            %OLDER PROCESSING FIX:
            if exist(obj.Params.T1toDWI.in.T1,'file') == 0
                [~ , obj.Params.T1toDWI.in.T1 ] = system(['ls ' obj.session_location 'T1' filesep 'gnc_T1.nii' ]);
                obj.Params.T1toDWI.in.T1=strtrim(obj.Params.T1toDWI.in.T1);
            end
            obj.proc_T1toDWI(); % ~~~> NOT IN COMPILED VERSION
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For DWI combined DTIFIT:
            %We will use the --wls option as it seems to improve the fit of
            %the diffusion tensor model and negativity values due to noise
            %REF: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=FSL;735f6320.1409
            %ALSO WE WILL ONLY USE THE 2.5k because w/ 7.5k the values
            %diverge for linearity (2.5k is close to linearity than
            %non-linear signal drop at 7.5k)
            obj.Params.Dtifit.in.movefiles = [ '..' filesep '06_2k5DTIFIT' ];
            obj.Params.Dtifit.in.fn ={obj.Params.CoRegMultiple.out.fn{1}}; %idx=1 for set4, lowest DWI gradient
            obj.Params.Dtifit.in.prefix = 'DTIFIT_FSLv509' ; %Double check this so you prefix the version of FSL!
            obj.Params.Dtifit.in.bvecs = {obj.Params.CoRegMultiple.out.bvecs{1}};
            obj.Params.Dtifit.in.bvals = {obj.Params.CoRegMultiple.out.bvals{1}};
            obj.Params.Dtifit.in.mask = {obj.Params.CoRegMultiple.out.combined_mask};
            
            obj.proc_dtifit();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %For GQI:
            obj.Params.GQI.in.movefiles = [ '..' filesep '07_Combined_Recon_gqi' ];
            obj.Params.GQI.in.fn = {obj.Params.CoRegMultiple.out.combined_fn};
            obj.Params.GQI.in.mask ={obj.Params.CoRegMultiple.out.combined_mask};
            obj.Params.GQI.in.bvecs = {obj.Params.CoRegMultiple.out.combined_bvecs};
            obj.Params.GQI.in.bvals = {obj.Params.CoRegMultiple.out.combined_bvals};
            
            obj.Params.GQI.in.prefix = 'GQI_DSISv053117' ; %Double check this so you prefix the version of DSISTUDIO!
            obj.Params.GQI.out.export = 'gfa,nqa0,nqa1';
            
            obj.Params.GQI.in.method = '4';    %for gqi model
            obj.Params.GQI.in.num_fiber = '3'; %modeling 3 fiber population
            obj.Params.GQI.in.param0 = '1.25'; %default parameter for gqi
            
            obj.proc_gqi();
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %FS2dwi:
            obj.Params.FS2dwi.in.movefiles = ['..' filesep '07_Combined_FS2dwi' ];
            obj.Params.FS2dwi.in.b0 =  {obj.Params.CoRegMultiple.out.combined_bet}; % Removed due to nonskulls stripping --> {obj.Params.CoRegMultiple.out.combined_b0} ;
            obj.Params.FS2dwi.in.aparcaseg = obj.Params.FreeSurfer.out.aparcaseg;
            obj.Params.FS2dwi.in.aparcaseg2009 = ...
                strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','aparc.a2009s+aseg'));
            
            %Using this will allow us to automatically select hippo-fields without discriminating whetehr they come from only a T1 or using T1-T2:            
            [~ , tmp_hippo_L ] = system(['ls ' strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','lh.hippoSfLabels*') ' | tail -1 ']);
            obj.Params.FS2dwi.in.hippofield_left = strtrim(tmp_hippo_L);
            [~ , tmp_hippo_R ] = system(['ls ' strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','rh.hippoSfLabels*') ' | tail -1 ']);
            obj.Params.FS2dwi.in.hippofield_right = strtrim(tmp_hippo_R);
            clear tmp_hippo_L tmp_hippo_R;
            %             obj.Params.FS2dwi.in.hippofield_left = ...
            %                 strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','lh.hippoSfLabels-T1-T2.v10'));
            %             obj.Params.FS2dwi.in.hippofield_right = ...
            %                 strtrim(strrep(obj.Params.FreeSurfer.out.aparcaseg,'aparc+aseg','rh.hippoSfLabels-T1-T2.v10'));
            obj.Params.FS2dwi.in.tmpfile_aparcaseg = [ obj.dependencies_dir filesep 'FS_DEPS' filesep  'FS_aparc.txt' ];
            obj.Params.FS2dwi.in.tmpfile_aparcaseg2009 = [ obj.dependencies_dir filesep 'FS_DEPS' filesep   'FS_aparc2009.txt' ];
            obj.Params.FS2dwi.in.tmpfile_hippo_bil = [ obj.dependencies_dir filesep 'FS_DEPS' filesep   'FS_hippolabels_bil.txt' ];
            
            obj.proc_FS2dwi();
            
            %~~~~> goint now to obj.CommonPostProc();
        end
        
        function obj = CommonPostProc(obj)
            %PREPARE T1 for Normalization using spm:
            obj.prep_spmT1_proc(); % ~~~~> NOT INCLUDED IN COMPILED
            %VERSION! 
         
            %TRACULA RELATED:
            for tohide=1:1
                obj.Params.Tracula.in.movefiles = ['..' filesep 'post_TRACULA' ];
                obj.Params.Tracula.in.fn = obj.Params.CoRegMultiple.out.combined_fn;
                obj.Params.Tracula.in.dcmrirc = [obj.dependencies_dir filesep 'TRACULA_DEPS' filesep 'dcmrirc.template' ];
                obj.Params.Tracula.in.FSDIR = obj.Params.FreeSurfer.dir;
                obj.Params.Tracula.in.bvec = obj.Params.CoRegMultiple.out.combined_bvecs;
                obj.Params.Tracula.in.bval = obj.Params.CoRegMultiple.out.combined_bvals;
                obj.Params.Tracula.in.nb0 = 28;
                
                obj.proc_tracula();
            end
            
            %TRKLAND RELATED - FX:
            for tohide=1:1
                obj.Trkland.root = [ obj.root  'post_TRKLAND' filesep ];
                %b0 params:
                obj.Trkland.fx.in.b0 = obj.Params.CoRegMultiple.out.combined_b0;
                obj.Trkland.fx.in.FA = obj.Params.Dtifit.out.FA{1};
                %Based on orientation, we will chooose a specific template (usually RAS)
                [~, Ori ] = system([obj.init_FS filesep 'bin' filesep 'mri_info ' obj.Trkland.fx.in.b0 ' | grep Orientation | awk ''{print $3}'''] );
                obj.Trkland.fx.tmp.ori = strtrim(Ori); clear Ori;
        
                %Fib params:
                obj.Trkland.fx.in.fib =strtrim(obj.Params.GQI.out.fibs_fn{end});
                if exist(obj.Trkland.fx.in.fib) == 0 ; error('No fib found in variable: trkland.trks.fx.in.fib. Please check!') ; end
                
                %Interpolation n:
                obj.Trkland.fx.in.n_interp=40; %According to ~average value on previous studies in connectome!
             
                obj.trkland_fx();  %~~~> NOT INCLUDED IN COMPILED
               % VERSION!!
            end
          
            %MODIFIED OR DEPRECATED CODE, CHECK COMENTED CODE BELOW:
            %% [ TO MODIFY OR DEPRECATED METHODS ]
       
            
            for tomodify_or_deprecated=1:1
                %% [HIDDEN]
                % Multi-shell bedpostX QBOOT 
                %          (not in used for now as we need at least 3 
                %          b-values to work correctly. Results were compared 
                %          and seemed to work but I decided to keep this
                %          option hidden, check: 
                %          https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/UserGuide#qboot_-_Estimation_of_fibre_orientations_using_q-ball_ODFs_and_residual_bootstrap 
                %          for more details)
                for tohide =1:1
                    %At this time and for this project this will be avoided since Qboot ideally works
                    %with three b-values instead of only two.
                    %For more information, check:
                    %https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/UserGuide#qboot_-_Estimation_of_fibre_orientations_using_q-ball_ODFs_and_residual_bootstrap
                    %             obj.Params.Qboot.in.movefiles = ['..' filesep 'post_Qboot' ];
                    %             obj.Params.Qboot.in.fn = obj.Params.CoRegMultiple.out.combined_fn;
                    %             obj.Params.Qboot.in.bvec = obj.Params.CoRegMultiple.out.combined_bvecs;
                    %             obj.Params.Qboot.in.bval = obj.Params.CoRegMultiple.out.combined_bvals;
                    %
                    %obj.proc_qboot();
                end
            
                %% [INACTIVE TRKLAND_CINGULUM TRKLAND_HIPPOCING]
                %TRKLAND yields in tractography hollowed Regions, at this
                %moment we only have done this for the fornix and previous
                %implementation of TRKLAND_HIPPOCING AND TRKLAND_CINGULUM
                %contains dilarted FreeSurfer parcellations for seed
                %tractography, not ideal as it does not follow the
                %TRKLAND_FX methods. For more info, check
                %rotrk_trimmedbyTOI
                % [INACTIVE ROIs/ROAs] TRKLAND_HIPPOCING:
                for tohide=1:1
                    %                 obj.Trkland.hippocing.in.hippo_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Left-Hippocampus.nii.gz');
                    %                 obj.Trkland.hippocing.in.hippo_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc2009_aseg/dwi_fs_Right-Hippocampus.nii.gz');
                    %                 obj.Trkland.hippocing.in.postcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-posteriorcingulate.nii.gz');
                    %                 obj.Trkland.hippocing.in.postcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-posteriorcingulate.nii.gz');
                    %
                    %                 %Interpolation n (for cingulum):
                    %                 obj.Trkland.hippocing.in.n_interp=33;
                    %                 obj.trkland_hippocing();
                %  [INACTIVE ROIs/ROAs] TRKLAND_CINGULUM:
                    %                 obj.Trkland.cingulum.in.rostantcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-rostralanteriorcingulate.nii.gz');
                    %                 obj.Trkland.cingulum.in.rostantcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-rostralanteriorcingulate.nii.gz');
                    %                 obj.Trkland.cingulum.in.postcing_lh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-lh-posteriorcingulate.nii.gz');
                    %                 obj.Trkland.cingulum.in.postcing_rh = strrep(obj.Params.FS2dwi.out.fn_aparc,'dwi_aparc+aseg.nii.gz','aparc_aseg/dwi_ctx-rh-posteriorcingulate.nii.gz');
                    %
                    %                 % BASE ON THESE VALUES (FROM EARLIER ADRC_PROCESSING W/O INTERP:
                    %                 % ninter_fx = 40;
                    %                 % ninter_cingulum = 32;
                    %                 % ninter_hippocing = 33;
                    %                 obj.Trkland.cingulum.in.n_interp = 32;
                    %                 trkland_cingulum(obj);
                end
                
                %% [ DEPRECATED METHODS ]
                for tohide_deprecated=1:1
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
                end
            end
        end
         
        %OTHER METHODS"
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
            obj.Params.tracxBYmask.allmasks.(tmp_txtfname).in.movefiles = ...
                ['..' filesep '..' filesep '..' filesep 'post_tracx' filesep 'all_masks' filesep  tmp_txtfname ];
            
            obj.Params.tracxBYmask.allmasks.(tmp_txtfname).probtracx2_args = ...
                ' -l --onewaycondition -c 0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --forcedir --opd  ' ;
            
            proc_tracxBYmask(obj,tmp_txtfname); %obj.resave()
            %%%%%%%%%%%%%%% END VARIABLE INITIALIZATION%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%% IMPLEMENTATION STARTS HERE %%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%% END OF IMPLEMENTATION  %%%%%%%%%%%%%%%%%%%%%
            
        end
        
        function obj = prep_spmT1_proc(obj)
            obj.Params.spmT1_Proc.in.tpm = [ fileparts(which('spm')) filesep 'tpm/TPM.nii'];
            if exist(obj.Params.spmT1_Proc.in.tpm,'file') == 0 
                error(['Cannot find: ' obj.Params.spmT1_Proc.in.tpm  ' Please check, exiting...' ]);
            end
            obj.proc_t1_spm();
        end
        
        function obj = proc_gradient_nonlin_correct(obj)
            wasRun = false;

            target = obj.Params.GradNonlinCorrect.in.target;
            [m h] = openIMG(target); if h.pinfo(1)~=1; h.pinfo(1)=1; spm_write_vol(h,m); end
            [a b c] = fileparts(target);
            outpath = obj.getPath(a,obj.Params.GradNonlinCorrect.in.movefiles);

            
            % addpath(genpath('/autofs/space/kant_002/users/rperea/DrigoScripts/adrc_diff/adrc_diff_prep/'));
            % mris_gradient_nonlin__unwarp_volume__batchmode_HCPS_v3(target, [outpath 'gc_mean.nii'], 'coeff_AS302.grad');
            
            infile = target;
            outfile = [outpath 'gnc_' b c];
            gradfile = obj.Params.GradNonlinCorrect.in.gradfile;
            
            
            %%% Compute the grdient nonlinearity correction
            if exist([outpath b '_deform_grad_rel.nii'],'file')==0
                cmd=['sh /autofs/space/kant_004/users/ConnectomeScanner/Scripts/adrc_diff_prep/run_mris_gradient_nonlin__unwarp_volume__batchmode_ADRC_v3.sh ' ...
                    '/usr/pubsw/common/matlab/8.5 ' ...
                    infile ' ' outfile ' ' gradfile ' '];
                system(cmd);
                wasRun = true;
            end
            obj.Params.GradNonlinCorrect.out.warpfile = [outpath b '_deform_grad_rel.nii'];
            
            %%% Apply the correction to the mean image.
            if exist(outfile,'file')==0
                cmd = ['applywarp -i ' infile ' -r ' infile ' -o ' outfile ' -w ' obj.Params.GradNonlinCorrect.out.warpfile ' --interp=spline'];
                runFS(cmd,pwd,3);
                system(['gunzip ' outpath '*.gz']);
                wasRun = true;
            end
            obj.Params.GradNonlinCorrect.out.meanimage = outfile;
            
            %%% Apply correction to full dataset
            fn = obj.Params.GradNonlinCorrect.in.fn;
            for ii = 1:numel(fn);
                infile = fn{ii};
                [a b c] = fileparts(infile);
                outpath = obj.getPath(a,obj.Params.GradNonlinCorrect.in.movefiles);
                outfile = [outpath 'gnc_' b c];
                if exist(outfile,'file')==0
                    cmd = ['applywarp -i ' infile ' -r ' infile ' -o ' outfile ' -w ' obj.Params.GradNonlinCorrect.out.warpfile ' --interp=spline'];
                    runFS(cmd,pwd,3);
                    system(['gunzip ' outpath '*.gz']);
                    wasRun = true;
                end
                obj.Params.GradNonlinCorrect.out.fn{ii,1} = outfile;
            end
            % obj.UpdateHist(obj.Params.GradNonlinCorrect,'proc_gradient_nonlin_correct',obj.Params.GradNonlinCorrect.out.warpfile,wasRun);
        end
    end
    methods ( Access = protected )
        function obj = getDCM2nii(obj,torun)
            %For proc_DCM2NII:
            obj.Params.DCM2NII.specific_vols=68;
            obj.Params.DCM2NII.seq_names={ 'ep2d_diff_7p5k_set1E60' 'ep2d_diff_7p5k_set2E60' ...
                'ep2d_diff_7p5k_set3E60' 'ep2d_diff_2p5k_set4E60' };
      
            for ii=1:4 % 4 sets of DWIs in this project!
                obj.Params.DCM2NII.in(ii).fsl2std_param = '-1 0 0 250.199 \n0 1 0 250.199 \n0 0 -1 0 \n0 0 0 1';
                obj.Params.DCM2NII.in(ii).prefix = obj.Params.DCM2NII.seq_names(ii);
                
                %DUE TO EARLIER UNPACKS AND NEWER ONE, WE WOULD NEED TO
                %SPLIT THE INPUT SELECTION IN TWO: 1) IF the SCANLOG files
                %EXISTS (EARLIER UNPACK) vs 2) IF THE NEWER UNPACK EXISTS
                obj.Params.DCM2NII.scanlog = [ obj.session_location filesep 'LogFiles' ...
                    filesep 'scan.log' ] ; % CASE 1)
                
                [tmp_bvecs_ok, tmp_bvecs] = ... 
                    system(([ 'ls ' obj.session_location filesep 'other' ...
                    filesep  '*' obj.Params.DCM2NII.seq_names{ii}  '*.bvecs' ]) );
                
                %OLDER UNPACK
                if exist(obj.Params.DCM2NII.scanlog,'file')
                    [ ~ , obj.Params.DCM2NII.in(ii).nvols ] = system([ 'cat ' ...
                        obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names{ii} ...
                        ' | tail -1 | awk ''{ print $7 }'' ' ]);
                      [ ~ , obj.Params.DCM2NII.in(ii).first_dcmfiles ]= system([ 'cat ' ...
                        obj.Params.DCM2NII.scanlog ' | grep ' obj.Params.DCM2NII.seq_names{ii} ...
                        ' | tail -1 | awk ''{ print $8 }'' ' ]);
                %NEWER UNPACK (Feb 2018):
                elseif exist(strtrim(tmp_bvecs),'file')
                    obj.Params.DCM2NII.newUnpack = true;
                    [ ~ , obj.Params.DCM2NII.in(ii).nvols ] = system(['cat ' ...
                        strtrim(tmp_bvecs) ' | wc -l ' ]);
                    obj.Params.DCM2NII.rawDiff = [ obj.session_location filesep 'other' filesep ] ;
                    obj.Params.DCM2NII.in(ii).first_dcmfiles = [] ; %No need to same first_dcmfile since its already been unpacked with *.bvecs
                else
                    fprintf(['DIFFUSION DICOMS NOT COMPLETE!\n 1.(newer unpack) Cannot find bvecs (newer unpack) in: '  strtrim(tmp_bvecs) ... 
                        '\n 2.(older unpack, before 2018) No scan.log found (older unpack) in: '  obj.Params.DCM2NII.scanlog ]);
                    
                    fprintf(['\nPLEASE CHECK THE LOGBOOK TO SEE IF DIFFUSION SEQUENCES WERE COLLECTED\n']);
                    fprintf('Throwing an error now....');
                    
                    error('Please verify that all the diffusion DCM files exist. Exiting now...');
                end
                obj.Params.DCM2NII.in(ii).nvols=str2num(obj.Params.DCM2NII.in(ii).nvols);
                obj.Params.DCM2NII.out(ii).location = [ obj.root 'Orig' filesep ];
                obj.Params.DCM2NII.out(ii).fn = [  obj.Params.DCM2NII.out(ii).location  cell2char_rdp(obj.Params.DCM2NII.seq_names(ii)) '.nii.gz' ];
            end
            obj.proc_dcm2nii();
        end
    end
end

