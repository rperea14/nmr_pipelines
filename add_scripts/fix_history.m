
AA=dir('./1*')
addpath('/cluster/brutha/MATLAB_Scripts/Utilities');
%AA = AA(27:140)


%*MAKE SURE YOU CHANGE YOUR dwiMRI_Session.m branch to redo_history!!
%BEFORE CONTINUE


tic
%BAD: ii=124 (AD93) 165 (AD142, no gnc_T1)
%BAD MISSING FS: 166 167 168 169 170 171 
%for ii= [  166 167 168 169 170 171 ] 

pths=MyPaths('hcp');

for ii=1:1 %numel(AA)
    
     SUBJID = AA(ii).name;
     
     fprintf(['\n\n\n\n\n\n\n\n IN ITERATION: ' num2str(ii) ' ID: ' SUBJID ]);
     obj_ADRC{ii} = load( [ pths.funcdir SUBJID '/DWIs/' SUBJID '.mat'] );
     AA=1;
     
     old_history{ii} = obj_ADRC{ii}.obj.history;
     for pp=1:numel(old_history{ii})
         if strcmp(old_history{ii}{pp}.lastRun(1:17),'proc_trkland_fx()')
             keep_trkland{ii} = old_history{ii}{pp};
         end
         if strcmp(old_history{ii}{pp}.lastRun(1:12),'proc_qboot()')
             keep_qboot{ii} = old_history{ii}{pp};
         end
     end
     AA=1;
     system(['mv '  pths.funcdir SUBJID '/DWIs/' SUBJID '.mat  '  pths.funcdir SUBJID '/DWIs/' SUBJID  '_BAK_03292018.mat' ]);
     
      obj_ADRC{ii}.obj=dwi_ADRC(SUBJID);
      
      obj_ADRC{ii}.obj.history{numel(obj_ADRC{ii}.obj.history)+1} =  keep_qboot{ii};
      obj_ADRC{ii}.obj.history{numel(obj_ADRC{ii}.obj.history)+1} =  keep_trkland{ii};
      
      obj_ADRC{ii}.obj.redo_history = false ; 
      obj_ADRC{ii}.obj.resave();
end
toc
timo=toc;
display(['Elapsed time is: ' num2str(timo/60) ' minutos']);
