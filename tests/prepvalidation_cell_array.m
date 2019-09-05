function [err,data,maxerr] = test(data,opts)

Parameters.par1 = {'a','b','c'};
Parameters.par2 = {'d','e'};

output = prepvalidation(Parameters);

err = ~isequal(size(output),[2*3 2]);
data = [];
maxerr = [];

end