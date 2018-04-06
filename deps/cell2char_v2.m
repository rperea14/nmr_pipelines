function X = cell2char(C)

for ii = 1:numel(C)
    if isnumeric(C{ii})
        C{ii} = num2str(C{ii});
    end 
end

X = char(C);