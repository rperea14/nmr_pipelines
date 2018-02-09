obj.setSPM12;
                    obj = replaceObjText(obj,{oldroot},{newroot});
%change line 48 - Talk to Aaron about directory output. 
%NOT SURE ABOUT THE PREFIXES HERE: ls /autofs/cluster/brutha/HABS_IIC/Sessions/170313_4PR00001/Structural/T1_gradFiles/ 

classdef dwi_HABSIIC < fMRI_Session
    
    properties
       root = '/autofs/cluster/brutha/HABS_IIC/Sessions/';
    end
    
    methods
        function obj = HABS_IIC_fcMRI(sessionname,opt)
            %%% Still not sure if the whole saving the objects as a mat
            %%% file thing is going to work.
            if nargin>1
                obj.root = opt;
            end
            
            obj.sessionname = sessionname;
            obj.root = [obj.root sessionname '/DWI/'];
            obj.objectHome = obj.root;
            
            newroot = obj.root;
            oldroot = obj.root;
            
            obj.setSPM12;
            obj.dosave = true;

            if exist([obj.objectHome filesep sessionname '.mat'],'file')>0 
                load([obj.objectHome filesep sessionname '.mat']);
                oldroot = obj.root;
                obj.wasLoaded = true;
            else
                obj.setMyParams; 
            end
            
            if nargin>1
                if ~strcmpi(oldroot,newroot)
                    obj = replaceObjText(obj,{oldroot},{newroot});
                    obj.resave;
                end
            end
        end
        
        function obj=setMyParams(obj)
            
            obj.rawfiles = dir_wfp([obj.root 'Orig/*.nii']);
            obj.fsdir='/autofs/eris/bang/ADRC/FreeSurfer6.0/';
            obj.fsubj=obj.sessionname;
            obj.vox = [2 2 2];
            obj.TR = 0.800;
            obj.bb = [-78 -112 -70; 78 76 90];
            obj.interporder = 5;
                
        end
        
        function obj = ProcessStructural(obj)
            rp = [fileparts(obj.root(1:end-1)) '/Structural/'];
            %%% Perform Gradient nonlin correction
            
            ff = char(dir_wfp([rp '0*T1*4e.nii']));
            [a b c] = fileparts(ff);
            
            if ~exist([a filesep 'gnc_' b c],'file')
                obj.Params.GradNonlinCorrect.in.movefiles = '';
                obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
                obj.Params.GradNonlinCorrect.in.gradfile = '/autofs/cluster/brutha/HABS_IIC/Scripts/coeff_AS82.grad';
                obj.Params.GradNonlinCorrect.in.fn = [];
                
                obj.Params.GradNonlinCorrect.out.warpfile = [];
                obj.Params.GradNonlinCorrect.out.meannii = [];
                obj.Params.GradNonlinCorrect.out.fn = [];
                
                obj.Params.GradNonlinCorrect.in.target = ff;
                obj.proc_gradient_nonlin_correct;
                
                mkdir([rp 'T1_gradFiles']);
                system(['mv ' rp '*grad*.nii ' rp 'T1_gradFiles/']);
                system(['mv ' rp '*jacobian*.nii ' rp 'T1_gradFiles/']);
                system(['mv ' rp '*voxel_volumes_warp*.nii ' rp 'T1_gradFiles/']);
            end
            
            ff = char(dir_wfp([rp '0*T2*vNav.nii']));
            [a b c] = fileparts(ff);
            
            if ~exist([a filesep 'gnc_' b c],'file')
                obj.Params.GradNonlinCorrect.in.movefiles = '';
                obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
                obj.Params.GradNonlinCorrect.in.gradfile = '/autofs/cluster/brutha/HABS_IIC/Scripts/coeff_AS82.grad';
                obj.Params.GradNonlinCorrect.in.fn = [];
                
                obj.Params.GradNonlinCorrect.out.warpfile = [];
                obj.Params.GradNonlinCorrect.out.meannii = [];
                obj.Params.GradNonlinCorrect.out.fn = [];
                
                obj.Params.GradNonlinCorrect.in.target = ff;
                obj.proc_gradient_nonlin_correct;
                
                mkdir([rp 'T2_gradFiles']);
                system(['mv ' rp '*grad*.nii ' rp 'T2_gradFiles/']);
                system(['mv ' rp '*jacobian*.nii ' rp 'T2_gradFiles/']);
                system(['mv ' rp '*voxel_volumes_warp*.nii ' rp 'T2_gradFiles/']);
            end
            %%%
            ff = dir_wfp([rp 'gnc*T1*.nii']);
            obj.Params.Coreg.in.source = ff{1};
            obj.Params.Coreg.in.target = '/autofs/cluster/brutha/MATLAB_Scripts/Atlas/HAB_T1_Avg.nii';
            obj.Params.Coreg.in.movecrfiles = [rp 'T1reg'];
            obj.Params.Coreg.in.movefiles = '';
            obj.Params.Coreg.in.style = 'spm';
            obj.Params.Coreg.in.getLabels = 0;
            obj.proc_coreg;
            %%%
            ff = dir_wfp([rp 'gnc*T2*.nii']);
            obj.Params.Coreg.in.source = ff{1};
            obj.Params.Coreg.in.target = obj.Params.Coreg.out.regimage;
            obj.Params.Coreg.in.movecrfiles = [rp 'T2reg'];
            obj.Params.Coreg.in.movefiles = '';
            obj.Params.Coreg.in.style = 'spm';
            obj.Params.Coreg.in.getLabels = 0;
            obj.proc_coreg;
            %%%
            
            if ~exist([rp 't1/t1.nii'],'file')
                mkdir([rp 't1']);
                ff = dir_wfp([rp 'cr_gnc_*T1*.nii']);
                copyfile(ff{1}, [rp 't1/t1.nii']);
            end
            
            obj.Params.spmT1_Proc.in.t1 = [rp 't1/t1.nii'];
            obj.Params.spmT1_Proc.in.outdir = [rp 't1/'];
            obj.proc_t1_spm;
            %%%
            
            if ~exist([rp 't2/t2.nii'],'file')
                mkdir([rp 't2']);
                ff = dir_wfp([rp 'cr_gnc_*T2*.nii']);
                copyfile(ff{1}, [rp 't2/t2.nii']);
            end
            
            obj.Params.spmT1_Proc.in.t1 = [rp 't2/t2.nii'];
            obj.Params.spmT1_Proc.in.outdir = [rp 't2/'];
            obj.proc_t1_spm;
            %%%
            if ~exist([rp 'combined/mwc6combo.nii'],'file');
                mkdir([rp 'combined']);
                
                prefix = {'c' 'wc' 'rc' 'mwc'};
                for zz = 1:numel(prefix)
                    
                    f1 = dir_wfp([rp 't1/' prefix{zz} '*t1.nii']);
                    f2 = dir_wfp([rp 't2/' prefix{zz} '*t2.nii']);
                    
                    for ii = 1:numel(f1);
                        [m1 h1] = openIMG(f1{ii});
                        [m2 h2] = openIMG(f2{ii});
                        
                        h = h1;
                        h.mat = (h1.mat+h2.mat)/2;
                        h.fname = [rp 'combined/' prefix{zz} num2str(ii) 'combo.nii'];
                        
                        spm_write_vol(h,(m1+m2)/2);
                    end
                end
            end
            %%%
            if ~exist([rp 'combined/iy_combo.nii'],'file');
                h1 = nifti([rp 't1/y_t1.nii']);
                h2 = nifti([rp 't2/y_t2.nii']);
                m = (h1.dat(:,:,:,:,:)+h2.dat(:,:,:,:,:))/2;
                
                h1.dat.fname = [rp 'combined/y_combo.nii'];
                create(h1);
                h1.dat(:,:,:,:,:) = m;
                clear h1
                %%%
                h1 = nifti([rp 't1/iy_t1.nii']);
                h2 = nifti([rp 't2/iy_t2.nii']);
                
                m = (h1.dat(:,:,:,:,:)+h2.dat(:,:,:,:,:))/2;
                
                h1.dat.fname = [rp 'combined/iy_combo.nii'];
                create(h1);
                h1.dat(:,:,:,:,:) = m;
                clear h1
                %%%
                
                [m1 h1] = openIMG([rp 't1/t1.nii']);
                [m2 h2] = openIMG([rp 't2/t2.nii']);
                newmat = (h1.mat+h2.mat)/2;
                h1.mat = newmat;
                h1.fname = [rp 'combined/t1.nii'];
                spm_write_vol(h1,m1);
                h1.fname = [rp 'combined/t2.nii'];
                spm_write_vol(h1,m2);
                
                [m1 h1] = openIMG([rp 't1/mt1.nii']);
                [m2 h2] = openIMG([rp 't2/mt2.nii']);
                newmat = (h1.mat+h2.mat)/2;
                h1.mat = newmat;
                h1.fname = [rp 'combined/mt1.nii'];
                spm_write_vol(h1,m1);
                h1.fname = [rp 'combined/mt2.nii'];
                spm_write_vol(h1,m2);
            end
            
            obj.Params.spmT1_Proc.in.outdir = [rp 'combined/'];
            obj.Params.spmT1_Proc.in.t1 = [rp 'combined/t1.nii'];
            obj.Params.spmT1_Proc.out.regfile = [rp 'combined/y_combo.nii'];
            obj.Params.spmT1_Proc.out.iregfile = [rp 'combined/iy_combo.nii'];
            obj.Params.spmT1_Proc.out.estTPM = dir_wfp([rp 'combined/c*.nii']);
            obj.resave;
        end
              
        function obj = CommonProc(obj)
            obj.ProcessStructural;
            
            rp = obj.root;
            
            if ~exist([rp 'Orig'],'dir')
                mkdir([rp 'Orig']);
                system(['mv ' rp '*.nii ' rp '/Orig/'])
            end
            %%%
            if isempty(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*REST_AP.nii']))
            
                obj.Params.GradNonlinCorrect.in.movefiles = '';
                obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
                obj.Params.GradNonlinCorrect.in.gradfile = '/autofs/cluster/brutha/HABS_IIC/Scripts/coeff_AS82.grad';
                fn = dir_wfp([rp 'Orig/0*_rfMRI_REST_AP.nii']);
                obj.Params.GradNonlinCorrect.in.fn = fn;
                
                mkdir([rp 'CommonProc/01_GNC/AP_epi_grads']);
                
                [m h] = openIMG(fn{1});
                h = h(1);
                h.fname = [rp 'CommonProc/01_GNC/AP_epi_grads/mean_AP.nii'];
                spm_write_vol(h,nanmean(m,4));
                
                obj.Params.GradNonlinCorrect.in.target = h.fname;

                obj.Params.GradNonlinCorrect.out.warpfile = [];
                obj.Params.GradNonlinCorrect.out.meannii = [];
                obj.Params.GradNonlinCorrect.out.fn = [];
                obj.proc_gradient_nonlin_correct;
                
                system(['mv ' rp 'Orig/gnc* ' rp 'CommonProc/01_GNC/']);
                
                %%%
                mkdir([rp 'CommonProc/01_GNC/AP_se_grads'])
                system(['cp ' rp 'Orig/0*SpinEcho*AP*.nii ' rp 'CommonProc/01_GNC/AP_se_grads/'])
                obj.Params.GradNonlinCorrect.in.fn = [];
                obj.Params.GradNonlinCorrect.in.target = char(dir_wfp([rp 'CommonProc/01_GNC/AP_se_grads/0*.nii']));
                obj.proc_gradient_nonlin_correct;
                
                system(['mv ' rp 'CommonProc/01_GNC/AP_se_grads/gnc_0*AP.nii ' rp 'CommonProc/01_GNC/']);
            end
            %%%
            if isempty(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*REST_PA.nii']))
            
                obj.Params.GradNonlinCorrect.in.movefiles = '';
                obj.Params.GradNonlinCorrect.in.prefix = 'gnc_';
                obj.Params.GradNonlinCorrect.in.gradfile = '/autofs/cluster/brutha/HABS_IIC/Scripts/coeff_AS82.grad';
                fn = dir_wfp([rp 'Orig/0*_rfMRI_REST_PA.nii']);
                obj.Params.GradNonlinCorrect.in.fn = fn;
                
                mkdir([rp 'CommonProc/01_GNC/PA_epi_grads']);
                
                [m h] = openIMG(fn{1});
                h = h(1);
                h.fname = [rp 'CommonProc/01_GNC/PA_epi_grads/mean_PA.nii'];
                spm_write_vol(h,nanmean(m,4));
                
                obj.Params.GradNonlinCorrect.in.target = h.fname;

                obj.Params.GradNonlinCorrect.out.warpfile = [];
                obj.Params.GradNonlinCorrect.out.meannii = [];
                obj.Params.GradNonlinCorrect.out.fn = [];
                obj.proc_gradient_nonlin_correct;
                
                system(['mv ' rp 'Orig/gnc* ' rp 'CommonProc/01_GNC/']);
                
                %%%
                mkdir([rp 'CommonProc/01_GNC/PA_se_grads'])
                system(['cp ' rp 'Orig/0*SpinEcho*PA*.nii ' rp 'CommonProc/01_GNC/PA_se_grads/'])
                obj.Params.GradNonlinCorrect.in.fn = [];
                obj.Params.GradNonlinCorrect.in.target = char(dir_wfp([rp 'CommonProc/01_GNC/PA_se_grads/0*.nii']));
                obj.proc_gradient_nonlin_correct;
                
                system(['mv ' rp 'CommonProc/01_GNC/PA_se_grads/gnc_0*PA.nii ' rp 'CommonProc/01_GNC/']);
            end
            
            %%%
            if isempty(dir_wfp([rp 'CommonProc/01_GNC/AP_epi_grads/cr_gnc_mean_AP.nii']))
                obj.Params.Coreg.in.source = char(dir_wfp([rp 'CommonProc/01_GNC/AP_epi_grads/gnc_mean_AP.nii']));
                obj.Params.Coreg.in.target = char(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*SpinEchoFieldMap_AP.nii']));
                obj.Params.Coreg.in.movecrfiles = [rp 'CommonProc/01_GNC/Reg1'];
                obj.proc_coreg;
                
                h = spm_vol(char(dir_wfp([rp 'CommonProc/01_GNC/AP_epi_grads/cr_gnc_mean_AP.nii'])));
                ff = char(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*_rfMRI_REST_AP.nii']));
                spm_get_space(ff, h.mat);
                delete([ff(1:end-3) 'mat']);
            end
            %%%
            if isempty(dir_wfp([rp 'CommonProc/01_GNC/PA_epi_grads/cr_gnc_mean_PA.nii']))
                obj.Params.Coreg.in.source = char(dir_wfp([rp 'CommonProc/01_GNC/PA_epi_grads/gnc_mean_PA.nii']));
                obj.Params.Coreg.in.target = char(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*SpinEchoFieldMap_PA.nii']));
                obj.Params.Coreg.in.movecrfiles = [rp 'CommonProc/01_GNC/Reg2'];
                obj.proc_coreg;
                
                h = spm_vol(char(dir_wfp([rp 'CommonProc/01_GNC/PA_epi_grads/cr_gnc_mean_PA.nii'])));
                ff = char(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*_rfMRI_REST_PA.nii']));
                spm_get_space(ff, h.mat);
                delete([ff(1:end-3) 'mat']);
            end
            %%%
            obj.Params.Realign.in.pars.rtm = 1;
            obj.Params.Realign.in.movefiles = [rp 'CommonProc/02_Realign/'];
            obj.proc_realign(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*_rfMRI_REST_AP.nii']));
            obj.proc_realign(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*_rfMRI_REST_PA.nii']));
            %%%
            
            %%% perform motion based unwarping
            uweflags.order = [12 12];
            uweflags.regorder = 1;
            uweflags.lambda = 100000;
            uweflags.jm = 0;
            uweflags.fot = [4 5];
            uweflags.sot = [];
            uweflags.fwhm = 4;
            uweflags.rem = 1;
            uweflags.noi = 5;
            uweflags.exp_round = 'Average';
            
            
            uwrflags.interp = 4;
            uwrflags.wrap = [0 0 0];
            uwrflags.mask = 1;
            uwrflags.which = 2;
            uwrflags.mean = 1;
            uwrflags.prefix = 'mu';
            uwrflags.udc = 1;
            
            if isempty(dir_wfp([rp 'CommonProc/02_Realign/mu_gnc_0*_rfMRI_REST_AP.nii']))
                h = spm_vol(char(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*_rfMRI_REST_AP.nii'])));
                uweflags.sfP = [];
                ds = spm_uw_estimate(h,uweflags);
                
                dsfile = spm_file(h(1).fname, 'suffix','_muw', 'ext','.mat');
                save(dsfile,'ds', '-v6');
                spm_uw_apply(ds,uwrflags);
                %%%
                fn = dir_wfp([rp 'CommonProc/01_GNC/mugnc_0*_rfMRI_REST_AP.nii']);
                [a b c] = fileparts(fn{1});
                system(['mv ' fn{1} ' ' rp 'CommonProc/02_Realign/mu_' b(3:end) c]);
                
                
                fn = dir_wfp([rp 'CommonProc/01_GNC/meanmugnc_0*_rfMRI_REST_AP.nii']);
                [a b c] = fileparts(fn{1});
                system(['mv ' fn{1} ' ' rp 'CommonProc/02_Realign/mean_mu_' b(7:end) c]);
            end
            %%%
            if isempty(dir_wfp([rp 'CommonProc/02_Realign/mu_gnc_0*_rfMRI_REST_PA.nii']))
                h = spm_vol(char(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*_rfMRI_REST_PA.nii'])));
                uweflags.sfP = [];
                ds = spm_uw_estimate(h,uweflags);
                
                dsfile = spm_file(h(1).fname, 'suffix','_muw', 'ext','.mat');
                save(dsfile,'ds', '-v6');
                spm_uw_apply(ds,uwrflags);
                
                %%%
                fn = dir_wfp([rp 'CommonProc/01_GNC/mugnc_0*_rfMRI_REST_PA.nii']);
                [a b c] = fileparts(fn{1});
                system(['mv ' fn{1} ' ' rp 'CommonProc/02_Realign/mu_' b(3:end) c]);
                
                
                fn = dir_wfp([rp 'CommonProc/01_GNC/meanmugnc_0*_rfMRI_REST_PA.nii']);
                [a b c] = fileparts(fn{1});
                system(['mv ' fn{1} ' ' rp 'CommonProc/02_Realign/mean_mu_' b(7:end) c]);
            end
            
            %%% resize to data to spin echos
            if isempty(dir_wfp([rp 'CommonProc/02_Realign/rs_mu_gnc_*AP.nii']))
                targ = spm_vol(char(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*SpinEcho*AP.nii'])));
                vol = spm_vol(char(dir_wfp([rp 'CommonProc/02_Realign/mean_mu_gnc_*AP.nii'])));
                m = resizeVol2(vol,targ,[3 0]);
                
                [a b c] = fileparts(vol.fname);
                writeIMG(targ,m,[a filesep 'rs_' b c],[16 0]);
                
                vol = spm_vol(char(dir_wfp([rp 'CommonProc/02_Realign/mu_gnc_*AP.nii'])));
                m = zeros([vol(1).dim numel(vol)]);
                for ii = 1:numel(vol)
                    m(:,:,:,ii) = resizeVol2(vol(ii),targ,[3 0]);
                end
                
                [a b c] = fileparts(vol(1).fname);
                writeIMG(targ,m,[a filesep 'rs_' b c],[16 0]);
            end
            %%%
            if isempty(dir_wfp([rp 'CommonProc/02_Realign/rs_mu_gnc_*PA.nii']))
                targ = spm_vol(char(dir_wfp([rp 'CommonProc/01_GNC/gnc_0*SpinEcho*PA.nii'])));
                vol = spm_vol(char(dir_wfp([rp 'CommonProc/02_Realign/mean_mu_gnc_*PA.nii'])));
                m = resizeVol2(vol,targ,[3 0]);
                
                [a b c] = fileparts(vol.fname);
                writeIMG(targ,m,[a filesep 'rs_' b c],[16 0]);
                
                vol = spm_vol(char(dir_wfp([rp 'CommonProc/02_Realign/mu_gnc_*PA.nii'])));
                m = zeros([vol(1).dim numel(vol)]);
                for ii = 1:numel(vol)
                    disp([ii numel(vol)]);
                    m(:,:,:,ii) = resizeVol2(vol(ii),targ,[3 0]);
                end
                
                [a b c] = fileparts(vol(1).fname);
                writeIMG(targ,m,[a filesep 'rs_' b c],[16 0]);
            end
            
            %%%     
            if isempty(dir_wfp([rp 'CommonProc/03_TopUp/warps/testOut1.nii']))
                
                mkdir([rp 'CommonProc/03_TopUp/warps/'])
                system(['cp ' rp 'CommonProc/02_Realign/rs_mean_mu_gnc_* ' rp 'CommonProc/03_TopUp/warps/']);
                
                [m1 h1] = openIMG(char(dir_wfp([rp 'CommonProc/03_TopUp/warps/rs_mean_mu_gnc_*_rfMRI_REST_PA.nii'])));
                [m2 h2] = openIMG(char(dir_wfp([rp 'CommonProc/03_TopUp/warps/rs_mean_mu_gnc_*_rfMRI_REST_AP.nii'])));
                
                h = h1;
                h = h(1);
                h.fname = [rp 'CommonProc/03_TopUp/warps/BothEchos.nii'];
                spm_write_vol(h,m1);
                h.n = [2 1];
                spm_write_vol(h,m2);
                
                %%%
                info = [0  1 0 0.104
                    0 -1 0 0.104];
                
                save([rp 'CommonProc/03_TopUp/warps/warp_info.txt'], 'info', '-ascii');
                
                %%%
                pth = [rp 'CommonProc/03_TopUp/warps/'];
                runFS(['topup --imain=' pth 'BothEchos.nii --datain=' pth 'warp_info.txt --config=b02b0.cnf --out=' pth 'results1 '],pwd,7);
                
                
                runFS(['applytopup --imain=' h2.fname ' --inindex=2 --method=jac --datain=' pth 'warp_info.txt --topup=' pth 'results1 --out=' pth 'testOut1'],pwd,7);
                runFS(['applytopup --imain=' h1.fname ' --inindex=1 --method=jac --datain=' pth 'warp_info.txt --topup=' pth 'results1 --out=' pth 'testOut2'],pwd,7);
                system(['gunzip ' pth '*.gz']);
            end
            
            %%%
            if isempty(dir_wfp([rp 'CommonProc/03_TopUp/su_*.nii']))
                pth = [rp 'CommonProc/03_TopUp/warps/'];
                
                f1 = dir_wfp([rp 'CommonProc/02_Realign/rs_mu_gnc*rfMRI_REST_AP.nii']);
                [a b c] = fileparts(f1{1}); nfn = [rp 'CommonProc/03_TopUp/su_' b c];
                runFS(['applytopup --imain=' f1{1} ' --inindex=2 --method=jac --datain=' pth 'warp_info.txt --topup=' pth 'results1 --out=' nfn],pwd,7);
                
                f1 = dir_wfp([rp 'CommonProc/02_Realign/rs_mu_gnc*rfMRI_REST_PA.nii']);
                [a b c] = fileparts(f1{1}); nfn = [rp 'CommonProc/03_TopUp/su_' b c];
                runFS(['applytopup --imain=' f1{1} ' --inindex=1 --method=jac --datain=' pth 'warp_info.txt --topup=' pth 'results1 --out=' nfn],pwd,7);
                
                system(['gunzip ' rp 'CommonProc/03_TopUp/*.gz']);
            end
            
            %%% run IUP
            obj.Params.Implicit_Unwarp.in.movefiles = [rp 'CommonProc/04_IUP/AP/'];
            obj.proc_implict_unwarping([rp 'CommonProc/03_TopUp/warps/testOut1.nii'],  dir_wfp([rp 'CommonProc/03_TopUp/su_rs_mu_gnc_*_rfMRI_REST_AP.nii']));
            
            %%%
            obj.Params.Implicit_Unwarp.in.movefiles = [rp 'CommonProc/04_IUP/PA/'];
            obj.proc_implict_unwarping([rp 'CommonProc/03_TopUp/warps/testOut2.nii'],  dir_wfp([rp 'CommonProc/03_TopUp/su_rs_mu_gnc*_rfMRI_REST_PA.nii']));
            
            %%%
            fn = dir_wfp2([rp 'CommonProc/04_IUP/*/uw_su_rs_mu_gnc*.nii']);
            obj.Params.DropVols.in.dropVols = 1:10;
            obj.Params.DropVols.in.movefiles = [rp '01_DropVols/'];
            obj.proc_drop_vols(fn);
            %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.ApplyNormNew.in.movefiles = [rp '02_Normed/'];
            obj.Params.ApplyNormNew.in.regfile = obj.Params.spmT1_Proc.out.regfile;
            obj.Params.ApplyNormNew.in.prefix = 'nn2_';
            
%             fn = dir_wfp2([rp 'CommonProc/04_IUP/*/uw_mean.nii']);
%             obj.Params.ApplyNormNew.in.fn = fn(1);
%             obj.Params.ApplyNormNew.in.fn = fn(2);
%             obj.proc_applynorm_new;
            
%             obj.Params.ApplyNormNew.out.normmean = obj.Params.ApplyNormNew.out.fn{1};
            
            obj.Params.ApplyNormNew.in.fn = obj.Params.DropVols.out.fn;
            obj.proc_applynorm_new;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%
            obj.Params.Smooth.in.kernel = [4 4 4];
            obj.Params.Smooth.in.prefix = 'ss4_';
            obj.Params.Smooth.in.movefiles = '../03_Smoothed/';
            obj.proc_smooth(obj.Params.ApplyNormNew.out.fn);
            
            obj.genTBRmaps(obj.Params.Smooth.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.resave;
        end
        
        function obj = BasicClean(obj)
            obj.Params.ComputePhysioRegs.in.movefiles = [obj.root '04_BasicClean/Regs/'];
            obj.Params.ComputePhysioRegs.in.type = 'indirecttpm';
            obj.Params.ComputePhysioRegs.in.whichparts = 2:3;
            obj.Params.ComputePhysioRegs.in.weighted = 1;
            obj.Params.ComputePhysioRegs.in.threshold = NaN;
            obj.Params.ComputePhysioRegs.in.prinComps = 1;
            obj.Params.ComputePhysioRegs.in.nPC = 10;
            obj.Params.ComputePhysioRegs.in.resample = [1 0];
            obj.Params.ComputePhysioRegs.in.masks = [];
            obj.Params.ComputePhysioRegs.in.filename = 'Mean_PhysioRegs2.txt';
            obj.proc_compute_physio_regs(obj.Params.DropVols.out.fn);
            % dir_wfp2([obj.root '/CommonProc/04_IUP/*/uw_su_rs_mu_gnc_*.nii'])
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.Realign.out.realigpars = dir_wfp([obj.root '/CommonProc/02_Realign/rp*.txt']);
            
            obj.Params.Filter.in.highcut = 0.01;
            obj.Params.Filter.in.lowcut = nan;
            
            obj.Params.CleanData.in.movefiles = [obj.root '/04_BasicClean/'];
            obj.Params.CleanData.in.filter = 1;
            obj.Params.CleanData.in.motion = 1;
            obj.Params.CleanData.in.physio = 1;
            obj.Params.CleanData.physioSquare = 0;
            obj.Params.CleanData.physioDeriv = 0;
            obj.Params.CleanData.motionSquare = 1;
            obj.Params.CleanData.motionDeriv = 1;
            obj.Params.CleanData.in.deriv  = 0;
            obj.Params.CleanData.in.square = 0;
            obj.Params.CleanData.in.reduce = 1;
            obj.proc_regress_clean(obj.Params.Smooth.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
            obj.genTBRmaps(obj.Params.CleanData.out.fn);
            obj.ExtractBOIs_v2(obj.Params.CleanData.out.fn);
            
            % tmp = table2array(D);
            % [rm pc s] = pca(zscore(tmp));
            % r = partialcorr(tmp,pc(:,1));
            % r = triu(r,1); r = r(r~=0);
            % hist(r,30); shg
        end
        
        function obj = LightClean(obj)
            obj.Params.Filter.in.highcut = [];
            obj.Params.Filter.in.lowcut = 0.01;
            
            obj.Params.CleanData.in.movefiles = '../LightClean/06_Cleaned/';
            obj.Params.CleanData.in.filter = 1;
            obj.Params.CleanData.in.motion = 1;
            obj.Params.CleanData.in.physio = 0;
            obj.Params.CleanData.in.deriv  = 1;
            obj.Params.CleanData.in.square = 1;
            obj.Params.CleanData.in.reduce = 0;
            obj.proc_regress_clean2(dir_wfp([obj.root '05_ImpUnwarp/uw_gnc_rr*.nii']));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.ApplyNormNew.in.fn = obj.Params.CleanData.out.fn;
            obj.Params.ApplyNormNew.in.movefiles = '../07_Normed/';
            obj.Params.ApplyNormNew.in.regfile = obj.Params.spmT1_Proc.out.regfile;
            obj.Params.ApplyNormNew.in.prefix = 'nn2_';
            obj.proc_applynorm_new;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.Smooth.in.kernel = [6 6 6];
            obj.Params.Smooth.in.prefix = 'ss6_';
            obj.Params.Smooth.in.movefiles = '../08_Smoothed/';
            obj.proc_smooth(obj.Params.ApplyNormNew.out.fn);
            
            %obj.genTBRmaps(obj.Params.Smooth.out.fn);
        end
        
        function obj = MediumClean(obj)
            obj.Params.ComputePhysioRegs.in.movefiles = 'PhysioRegs/';
            obj.Params.ComputePhysioRegs.in.type = 'indirecttpm';
            obj.Params.ComputePhysioRegs.in.whichparts = 1:3;
            obj.Params.ComputePhysioRegs.in.weighted = 1;
            obj.Params.ComputePhysioRegs.in.threshold = NaN;
            obj.Params.ComputePhysioRegs.in.prinComps = 0;
            obj.Params.ComputePhysioRegs.in.nPC = NaN;
            obj.Params.ComputePhysioRegs.in.resample = [1 0];
            obj.Params.ComputePhysioRegs.in.masks = [];
            obj.Params.ComputePhysioRegs.in.filename = 'Mean_PhysioRegs3.txt';
            obj.proc_compute_physio_regs(obj.Params.Implicit_Unwarp.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.Filter.in.highcut = 0;
            obj.Params.Filter.in.lowcut = 0.01;
            
            obj.Params.CleanData.in.movefiles = '../MediumClean/06_Cleaned/';
            obj.Params.CleanData.in.filter = 1;
            obj.Params.CleanData.in.motion = 1;
            obj.Params.CleanData.in.physio = 1;
            obj.Params.CleanData.in.deriv  = 1;
            obj.Params.CleanData.in.square = 1;
            obj.Params.CleanData.in.reduce = 0;
            obj.proc_regress_clean2(obj.Params.Implicit_Unwarp.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.ApplyNormNew.in.fn = obj.Params.CleanData.out.fn;
            obj.Params.ApplyNormNew.in.movefiles = '../07_Normed/';
            obj.Params.ApplyNormNew.in.regfile = obj.Params.spmT1_Proc.out.regfile;
            obj.Params.ApplyNormNew.in.prefix = 'nn2_';
            obj.proc_applynorm_new;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.Smooth.in.kernel = [6 6 6];
            obj.Params.Smooth.in.prefix = 'ss6_';
            obj.Params.Smooth.in.movefiles = '../08_Smoothed/';
            obj.proc_smooth(obj.Params.ApplyNormNew.out.fn);
            
            %obj.genTBRmaps(obj.Params.Smooth.out.fn);
        end
        
        function obj = MediumClean2(obj)
            obj.Params.ComputePhysioRegs.in.movefiles = 'PhysioRegs/';
            obj.Params.ComputePhysioRegs.in.type = 'indirecttpm';
            obj.Params.ComputePhysioRegs.in.whichparts = 2:3;
            obj.Params.ComputePhysioRegs.in.weighted = 1;
            obj.Params.ComputePhysioRegs.in.threshold = NaN;
            obj.Params.ComputePhysioRegs.in.prinComps = 1;
            obj.Params.ComputePhysioRegs.in.nPC = 10;
            obj.Params.ComputePhysioRegs.in.resample = [1 0];
            obj.Params.ComputePhysioRegs.in.masks = [];
            obj.Params.ComputePhysioRegs.in.filename = 'Mean_PhysioRegs2.txt';
            obj.proc_compute_physio_regs(obj.Params.Implicit_Unwarp.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.Filter.in.highcut = 0.01;
            obj.Params.Filter.in.lowcut = nan;
            
            obj.Params.CleanData.in.movefiles = '../MediumClean2/06_Cleaned/';
            obj.Params.CleanData.in.filter = 1;
            obj.Params.CleanData.in.motion = 1;
            obj.Params.CleanData.in.physio = 1;
            obj.Params.CleanData.in.deriv  = 1;
            obj.Params.CleanData.in.square = 1;
            obj.Params.CleanData.in.reduce = 1;
            obj.proc_regress_clean(obj.Params.Implicit_Unwarp.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.ApplyNormNew.in.fn = obj.Params.CleanData.out.fn;
            obj.Params.ApplyNormNew.in.movefiles = '../07_Normed/';
            obj.Params.ApplyNormNew.in.regfile = obj.Params.spmT1_Proc.out.regfile;
            obj.Params.ApplyNormNew.in.prefix = 'nn2_';
            obj.proc_applynorm_new;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.Smooth.in.kernel = [6 6 6];
            obj.Params.Smooth.in.prefix = 'ss6_';
            obj.Params.Smooth.in.movefiles = '../08_Smoothed/';
            obj.proc_smooth(obj.Params.ApplyNormNew.out.fn);
            
%             obj.Params.Smooth.in.kernel = [4 4 4];
%             obj.Params.Smooth.in.prefix = 'ss4_';
%             obj.Params.Smooth.in.movefiles = '../08_Smoothed/';
%             obj.proc_smooth(obj.Params.ApplyNormNew.out.fn);
            
            obj.genTBRmaps(obj.Params.Smooth.out.fn);
            obj.ExtractBOIs_v2(obj.Params.Smooth.out.fn)
        end
        
        function obj = TradClean(obj)
            obj.Params.ComputePhysioRegs.in.movefiles = 'PhysioRegs/';
            obj.Params.ComputePhysioRegs.in.type = 'indirecttpm';
            obj.Params.ComputePhysioRegs.in.whichparts = 1:3;
            obj.Params.ComputePhysioRegs.in.weighted = 1;
            obj.Params.ComputePhysioRegs.in.threshold = NaN;
            obj.Params.ComputePhysioRegs.in.prinComps = 0;
            obj.Params.ComputePhysioRegs.in.nPC = 10;
            obj.Params.ComputePhysioRegs.in.resample = [1 0];
            obj.Params.ComputePhysioRegs.in.masks = [];
            obj.Params.ComputePhysioRegs.in.filename = 'TradRegs.txt';
            obj.proc_compute_physio_regs(obj.Params.Implicit_Unwarp.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.Filter.in.highcut = 0.10;
            obj.Params.Filter.in.lowcut  = 0.01;
            
            obj.Params.CleanData.in.movefiles = '../TradClean/06_Cleaned/';
            obj.Params.CleanData.in.filter = 1;
            obj.Params.CleanData.in.motion = 1;
            obj.Params.CleanData.in.physio = 1;
            obj.Params.CleanData.in.deriv  = 1;
            obj.Params.CleanData.in.square = 0;
            obj.Params.CleanData.in.reduce = 0;
            obj.proc_regress_clean(obj.Params.Implicit_Unwarp.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.ApplyNormNew.in.fn = obj.Params.CleanData.out.fn;
            obj.Params.ApplyNormNew.in.movefiles = '../07_Normed/';
            obj.Params.ApplyNormNew.in.regfile = obj.Params.spmT1_Proc.out.regfile;
            obj.Params.ApplyNormNew.in.prefix = 'nn2_';
            obj.proc_applynorm_new;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.Smooth.in.kernel = [6 6 6];
            obj.Params.Smooth.in.prefix = 'ss6_';
            obj.Params.Smooth.in.movefiles = '../08_Smoothed/';
            obj.proc_smooth(obj.Params.ApplyNormNew.out.fn);
            
%             obj.Params.Smooth.in.kernel = [4 4 4];
%             obj.Params.Smooth.in.prefix = 'ss4_';
%             obj.Params.Smooth.in.movefiles = '../08_Smoothed/';
%             obj.proc_smooth(obj.Params.ApplyNormNew.out.fn);
            
            obj.genTBRmaps(obj.Params.Smooth.out.fn);
            obj.ExtractBOIs_v2(obj.Params.Smooth.out.fn)
        end
        
        function obj = BucknerStyle(obj)
            obj.Params.ApplyNormNew.in.fn = obj.Params.Implicit_Unwarp.out.fn;
            obj.Params.ApplyNormNew.in.movefiles = '../06_Normed/';
            obj.Params.ApplyNormNew.in.regfile = obj.Params.spmT1_Proc.out.regfile;
            obj.Params.ApplyNormNew.in.prefix = 'nn2_';
            obj.proc_applynorm_new;
            
            
            obj.Params.Smooth.in.kernel = [6 6 6];
            obj.Params.Smooth.in.prefix = 'ss6_';
            obj.Params.Smooth.in.movefiles = '../07_Smoothed/';
            obj.proc_smooth(obj.Params.ApplyNormNew.out.fn);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.Params.ComputePhysioRegs.in.movefiles = 'PhysioRegs/';
            obj.Params.ComputePhysioRegs.in.type = 'masks';
            obj.Params.ComputePhysioRegs.in.whichparts = [];
            obj.Params.ComputePhysioRegs.in.weighted = 0;
            obj.Params.ComputePhysioRegs.in.threshold = .1;
            obj.Params.ComputePhysioRegs.in.prinComps = false;
            obj.Params.ComputePhysioRegs.in.nPC = NaN;
            obj.Params.ComputePhysioRegs.in.resample = [0 0];
            obj.Params.ComputePhysioRegs.in.masks = dir_wfp('/autofs/cluster/brutha/MATLAB_Scripts/RestingState/avg152T1*.nii');
            obj.Params.ComputePhysioRegs.in.filename = 'BucknerRegs.txt';
            obj.proc_compute_physio_regs(obj.Params.Smooth.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            obj.Params.Filter.in.highcut = 0.08;
            obj.Params.Filter.in.lowcut  = 0.01;
            
            obj.Params.CleanData.in.movefiles = '../08_Cleaned/';
            obj.Params.CleanData.in.filter = 1;
            obj.Params.CleanData.in.motion = 1;
            obj.Params.CleanData.in.physio = 1;
            obj.Params.CleanData.in.deriv  = 1;
            obj.Params.CleanData.in.square = 0;
            obj.Params.CleanData.in.reduce = 0;
            obj.proc_regress_clean(obj.Params.Smooth.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.genTBRmaps(obj.Params.CleanData.out.fn);
            obj.ExtractBOIs_v2(obj.Params.CleanData.out.fn);
         end
        
        function obj = HeavyClean(obj)
            obj.Params.ComputePhysioRegs.in.movefiles = 'PhysioRegs/';
            obj.Params.ComputePhysioRegs.in.type = 'indirecttpm';
            obj.Params.ComputePhysioRegs.in.whichparts = 1:6;
            obj.Params.ComputePhysioRegs.in.weighted = 1;
            obj.Params.ComputePhysioRegs.in.threshold = NaN;
            obj.Params.ComputePhysioRegs.in.prinComps = 1;
            obj.Params.ComputePhysioRegs.in.nPC = 10;
            obj.Params.ComputePhysioRegs.in.resample = [1 0];
            obj.Params.ComputePhysioRegs.in.masks = [];
            obj.Params.ComputePhysioRegs.in.filename = 'PCA_PhysioRegs.txt';
            obj.proc_compute_physio_regs(obj.Params.Implicit_Unwarp.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.Filter.in.highcut = 0;
            obj.Params.Filter.in.lowcut = 0.01;
            
            obj.Params.CleanData.in.movefiles = '../HeavyClean/06_Cleaned/';
            obj.Params.CleanData.in.filter = 1;
            obj.Params.CleanData.in.motion = 1;
            obj.Params.CleanData.in.physio = 1;
            obj.Params.CleanData.in.deriv  = 1;
            obj.Params.CleanData.in.square = 1;
            obj.Params.CleanData.in.reduce = 1;
            obj.proc_regress_clean2(obj.Params.Implicit_Unwarp.out.fn);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.ApplyNormNew.in.fn = obj.Params.CleanData.out.fn;
            obj.Params.ApplyNormNew.in.movefiles = '../07_Normed/';
            obj.Params.ApplyNormNew.in.regfile = obj.Params.spmT1_Proc.out.regfile;
            obj.Params.ApplyNormNew.in.prefix = 'nn2_';
            obj.proc_applynorm_new;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.Params.Smooth.in.kernel = [6 6 6];
            obj.Params.Smooth.in.prefix = 'ss6_';
            obj.Params.Smooth.in.movefiles = '../08_Smoothed/';
            obj.proc_smooth(obj.Params.ApplyNormNew.out.fn);
            
%             obj.genTBRmaps(obj.Params.Smooth.out.fn);
        end
        
        function obj = FixClean(obj)
            
        end
        
        function obj = MapBOIsToNativeSpace(obj)
            obj.setSPM12;
            fn = dir_wfp('/autofs/space/aristotle_003/users/FC_and_Cognition/TBR_Templates/FROIS/renamed/*.nii');
            fn = [fn; dir_wfp('/autofs/space/aristotle_003/users/FC_and_Cognition/TBR_Templates/SensoryBOIs2/*.nii')];
            
            obj.Params.ApplyReverseNormNew.in.movefiles = [obj.root 'BOIs/'];
            obj.Params.ApplyReverseNormNew.in.fn = fn;
            obj.Params.ApplyReverseNormNew.in.targ = obj.Params.Implicit_Unwarp.out.newmean;
            obj.Params.ApplyReverseNormNew.in.regfile = obj.Params.spmT1_Proc.out.iregfile;
            
            obj.proc_apply_reservsenorm_new;
        end
        
        function obj = ExtractBOIs_v1(obj)
            
            [msk name] = dir_wfp([obj.root 'BOIs/*.nii']);
            fn = obj.Params.Implicit_Unwarp.out.fn;
            
            fields = regexprep(name,'in_|.nii','');
            obj.Params.BOI.V1.names = regexprep(regexprep(name,'in_|.nii',''),'_',' ');
            
            
            for ii = 1:numel(fn);
                disp(ii);
                dat = table;
                m = FastRead(fn{ii});
                
                for jj = 1:numel(msk)
                    dat.(fields{jj}) = nanmean(m(:,FastRead(msk{jj})==1),2);
                end
                
                obj.Params.BOI.V1.(['Run' sprintf('%0.2i',ii)]) = dat;
            end
            obj.resave;
        end
        
        function obj = ExtractBOIs_v2(obj,fn)
            [msk name] = dir_wfp('/autofs/cluster/brutha/MATLAB_Scripts/Atlas/fMRI/SchultzMaps/FROI/*.nii');
            
            fields = regexprep(name,'in_|.nii','');
            %names = regexprep(regexprep(name,'in_|.nii',''),'_',' ');
            
            D = [];
            h = spm_vol([fn{1} ',1']);
            for ii = 1:numel(fn);
                disp(ii);
                dat = table;
                m = FastRead(fn{ii});
                
                for jj = 1:numel(msk)
                    mm = resizeVol2(spm_vol(msk{jj}),h,[0 0]);
                    dat.(fields{jj}) = zscore(nanmean(m(:,mm==1),2));
                end
                D = [D; dat];
            end
            
            [a b c] = fileparts(fn{1});
            save([a '/FROI.mat'],'D');
        end
        
        function obj = genTBRmaps(obj,fn)
            tt = dir_wfp('/autofs/cluster/brutha/MATLAB_Scripts/Atlas/fMRI/SchultzMaps/StandardTemplates/*.nii');
            if ~exist([fileparts(fn{1}) filesep 'TBR_Maps_Standard'],'dir')
                tt = dir_wfp('/autofs/cluster/brutha/MATLAB_Scripts/Atlas/fMRI/SchultzMaps/StandardTemplates/*.nii');
                TBR(fn,tt,[],[],'_Standard',[fileparts(fn{1}) filesep]);
            end
            
            for ii = 1:numel(fn)
                if ~exist([fileparts(fn{ii}) filesep 'TBR_Maps_Standard_Run' sprintf('%0.2i',ii)],'dir')
                    TBR(fn{ii},tt,[],[],['_Standard_Run' sprintf('%0.2i',ii)],[fileparts(fn{ii}) filesep]);
                end
            end
            
            %%% Save paths to output files somewhere?
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
            
            obj.UpdateHist(obj.Params.GradNonlinCorrect,'proc_gradient_nonlin_correct',obj.Params.GradNonlinCorrect.out.warpfile,wasRun);
        end
        
        function resave(obj)
            save([obj.objectHome filesep obj.sessionname '.mat'],'obj');
        end
        
        function setSPM12(obj)
            if isempty(strfind(which('spm.m'),'spm12'))
                tmp = path;
                tmp = regexp(tmp,':','split')';
                tmp = tmp(searchCellStr('(/spm8/|spm8$)',tmp));
                
                for ii = 1:numel(tmp)
                    rmpath(tmp{ii});
                end
                
                addpath('/autofs/cluster/brutha/spm12/')
            end
        end
        
        function setSPM8(obj)
            if isempty(strfind(which('spm.m'),'spm9'))
                tmp = path;
                tmp = regexp(tmp,':','split')';
                tmp = tmp(searchCellStr('(/spm12/|spm12$)',tmp));
                
                for ii = 1:numel(tmp)
                    rmpath(tmp{ii});
                end
                
                addpath('/autofs/cluster/brutha/MATLAB_Scripts/spm8/')
            end
        end
    end
end