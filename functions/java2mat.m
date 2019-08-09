function varargout = java2mat(javaobject)

out = cell(javaobject);
varargout = cell(length(out),1);
for i=1:length(out)
   
    if isa(out{i},'java.lang.Object')
        out{i} = cell(out{i});
    end
    if isa(out{i},'numeric') && iscolumn(out{i})
        out{i} = out{i}.';
    end
    varargout{i} = out{i};
end

end