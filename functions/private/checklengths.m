function boolean = checklengths(arrayA,arrayB)

if length(arrayA) ~= length(arrayB)
    error('The ''%s'' and ''%s'' variables must be equally long.',inputname(1),inputname(2));
end
boolean = true;
end