
AA=dir('./1*');
addpath('/cluster/brutha/MATLAB_Scripts/Utilities');
%AA = AA(27:140)


%*MAKE SURE YOU CHANGE YOUR dwiMRI_Session.m branch to redo_history!!
%BEFORE CONTINUE


%** THIS WILL WORK ONLY In PRE_PROCESS STEPS! 

%Edited by Rodrigo Perea grandrigo@gmail.com 

tic
%BAD: ii=124 (AD93) 165 (AD142, no gnc_T1)
%BAD MISSING FS: 166 167 168 169 170 171 
%for ii= [  166 167 168 169 170 171 ] 

pths=MyPaths('hcp');

for ii=1:numel(AA)
    
     SUBJID = AA(ii).name;
     
     fprintf(['\n\n\n\n\n\n\n\n IN ITERATION: ' num2str(ii) ' ID: ' SUBJID ]);
     obj_ADRC{ii} = load( [ pths.funcdir SUBJID '/DWIs/' SUBJID '.mat'] );
     AA=1;
     
     old_history{ii} = obj_ADRC{ii}.obj.history;
    % keep_trkland{ii}='';
     keep_qboot{ii} = '' ;
     keep_T1toDWI{ii} = '' ;
     keep_t1_spm{ii} = '' ;
     
     for pp=1:numel(old_history{ii})
%          if strcmp(old_history{ii}{pp}.lastRun(1:17),'proc_trkland_fx()')
%              keep_trkland{ii} = old_history{ii}{pp};
%          end
         if strcmp(old_history{ii}{pp}.lastRun(1:12),'proc_qboot()')
             keep_qboot{ii} = old_history{ii}{pp};
         end
         
         if strcmp(old_history{ii}{pp}.lastRun(1:14),'proc_T1toDWI()')
             keep_T1toDWI{ii}= old_history{ii}{pp};
         end
         
          if strcmp(old_history{ii}{pp}.lastRun(1:13),'proc_t1_spm()')
             keep_t1_spm{ii}= old_history{ii}{pp};
          end
     end
     AA=1;
     system(['mv '  pths.funcdir SUBJID '/DWIs/' SUBJID '.mat  '  pths.funcdir SUBJID '/DWIs/' SUBJID  '_BAK_03292018.mat' ]);
     
      obj_ADRC{ii}.obj=dwi_ADRC(SUBJID);
      
      %Re-reading all the trks in TRKLAND
      
      %Keeping these histories:
      if ~isempty(keep_qboot{ii}) ; obj_ADRC{ii}.obj.history{numel(obj_ADRC{ii}.obj.history)+1} =  keep_qboot{ii}; end
      obj_ADRC{ii}.obj.trkland_fx();
      %if ~isempty(keep_trkland{ii}) ; obj_ADRC{ii}.obj.history{numel(obj_ADRC{ii}.obj.history)+1} =  keep_trkland{ii}; end
      if ~isempty(keep_T1toDWI{ii}) ; obj_ADRC{ii}.obj.history{numel(obj_ADRC{ii}.obj.history)+1} =  keep_T1toDWI{ii}; end
      if ~isempty(keep_t1_spm{ii}) ; obj_ADRC{ii}.obj.history{numel(obj_ADRC{ii}.obj.history)+1} =  keep_t1_spm{ii}; end
      
      obj_ADRC{ii}.obj.redo_history = false ; 
      obj_ADRC{ii}.obj.resave();
end
toc
timo=toc;
display(['Elapsed time is: ' num2str(timo/60) ' minutos']);
