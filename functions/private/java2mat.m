function varargout = java2mat(javaobject)

%Get the outputs cached on the Java hashmap
out = cell(javaobject);
%Prepare MATLAB output container
varargout = cell(length(out),1);

%Loop over chached outputs
for i=1:size(out,2)
    
    %If the output is numeric
    if isa(out{1,i},'numeric')
        %Reshape to appropiate chached shape
        shape = out{2,i}.';
        out{i} = reshape(out{1,i},shape);
    end
    %If the output is a linked hashmap reconvert to struct
    if isa(out{1,i},'java.util.LinkedHashMap')
        table = out{1,i};
        fieldnames = table.keySet;
        fieldnames = cell(fieldnames.toArray);
        structure = struct(); 
        for j=1:length(fieldnames)
            structure = setfield(structure,fieldnames{j},table.get(fieldnames{j}));
        end
        out{1} = structure;
    end
    %If output is a specific type of object, reconvert
    if isa(out{i},'java.lang.Object')
        out{1,i} = cell(out{1,i});
    end
    %Otherwise just store the output directly
    varargout{i} = out{1,i};
end

end