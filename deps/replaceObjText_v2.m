function object = replaceObjText_v2(object,old,new)

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


function [AllFields, FieldType] = getFieldList(object)
AllFields = [];
FieldType = [];
newfields = {'object'};
while true
    newfields2 = [];
    for ii = 1:numel(newfields)
        [thesefields, otherfields, types] = getFields2(eval(newfields{ii}),newfields{ii});
        
        AllFields = [AllFields; thesefields];
        newfields2 = [newfields2; otherfields];
        FieldType = [FieldType; types];
    end
    
    
    if isempty(newfields2)
        break
    else
        newfields = newfields2;
    end
end

AllFields = regexprep(AllFields,'^object.','');
end



function [thesefields, otherfields, types] = getFields2(S,start)

thesefields = [];
otherfields = [];
types = [];

if isstruct(S) && numel(S)>1
    for ii = 1:numel(S)
        otherfields{end+1,1} = [start '(' num2str(ii) ')'];
    end
    return
end

if iscell(S) && numel(S)>1
    for ii = 1:numel(S)
        otherfields{end+1,1} = [start '{' num2str(ii) '}'];
    end
    return
end

tmp = fields(S);

if ~isempty(start); start = [start '.']; end

for ii = 1:numel(tmp)
    try
        a = fields(S.(tmp{ii}));
        otherfields{end+1,1} = [start tmp{ii}];
    catch
        try; class([S.(tmp{ii})]); catch; continue; end;
        if strcmp(class([S.(tmp{ii})]),'cell')
            if  ~isempty(S.(tmp{ii})) && isstruct(S.(tmp{ii}){1})
                
                for jj = 1:numel(S.(tmp{ii}))
                    otherfields{end+1,1} = [start tmp{ii} '{' num2str(jj) '}'];
                end
            else
                thesefields{end+1,1} = [start tmp{ii}];
                try
                    x = char(S.(tmp{ii}));
                    types{end+1,1} = 'char';
                catch
                    types{end+1,1} = class(S.(tmp{ii}));
                end
            end
        else
            types{end+1,1} = class([S.(tmp{ii})]);
            thesefields{end+1,1} = [start tmp{ii}];
        end
    end
end
end