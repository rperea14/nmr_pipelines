AA=dir('./1*')
addpath('/cluster/brutha/MATLAB_Scripts/Utilities');
%AA = AA(27:140)

%(now ii shuold be 93)
tic
%BAD: ii=124 (AD93) 165 (AD142, no gnc_T1)
%BAD MISSING FS: 166 167 168 169 170 171 
%for ii= [  166 167 168 169 170 171 ] 

pths=MyPaths('hcp');

ref_obj=load('170718_8CSAD00151/DWIs/170718_8CSAD00151.mat');
replaced_obj = replaceObjText_v2(ref_obj.obj,{'170718_8CSAD00151'},{'140924_8CS00086'});
for ii=1:1 %numel(AA)
    
     SUBJID = AA(ii).name;

     
     
     
     fprintf(['\n\n\n\n\n\n\n\n IN ITERATION: ' num2str(ii) ' ID: ' SUBJID ]);
     obj_ADRC{ii} = load( [ pths.funcdir SUBJID '/DWIs/' SUBJID '.mat'] );
     
     
     
end
toc
timo=toc;
display(['Elapsed time is: ' num2str(timo/60) ' minutos']);




function object = rObjText_v2(object,old,new)

[AllFields, FieldType] = getFieldList(object);

i1 = searchCellStr('^char$',FieldType);
for ii =i1
    for jj = 1:numel(old)
        x = eval(['object.' AllFields{ii}]);
        ss=size(x);
        
        if ss(1) > 1 ; continue;end;
        if isempty(x); continue; end;
        if islogical(x); continue; end;
        if isnumeric(x); continue; end;
        
        
        %TODEBUG-->
        %display([num2str(ii) ' with size: ' num2str(ss(1)) ' ' num2str(ss(2)) ]);
        
        
        x1 = regexprep(x,old{jj},new{jj});
        
        if ~isequal(x,x1)
            eval(['object.' AllFields{ii} '=x1;']);
        end
    end
end
end

