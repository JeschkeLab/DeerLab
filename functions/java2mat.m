function varargout = java2mat(javaobject)

out = cell(javaobject);
varargout = cell(length(out),1);
for i=1:length(out)
    

    if isa(out{i},'numeric') && isvector(out{i}) && ~iscolumn(out{i}) && size()
        out{i} = out{i}.';
    end
    if isa(out{i},'java.util.Hashtable')
        table = out{i};
        fieldnames = table.keySet;
        fieldnames = cell(fieldnames.toArray);
        structure = struct(); 
        for j=1:length(fieldnames)
            structure = setfield(structure,fieldnames{j},table.get(fieldnames{j}));
        end
        out{1} = structure;
    end
    if isa(out{i},'java.lang.Object')
        out{i} = cell(out{i});
    end
    varargout{i} = out{i};
end

end