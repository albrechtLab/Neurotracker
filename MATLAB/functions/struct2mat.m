%% struct2mat 
%   pulls structure information into matrices

function output = struct2mat(dim,structure,index,fields)

    if ~iscell(fields) fields = {fields}; end
    if isempty(index) index = 1:numel(structure); end

    output = [];
    for i = 1:length(index)
        substructure = structure(index(i));
        for j = 1:length(fields)
            substructure = getfield(substructure,char(fields(j)));
        end

        output = safecat(dim,output,substructure);
    end
end