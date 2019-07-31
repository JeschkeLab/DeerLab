function [err,data,maxerr] = test(data,opts)

Parameters(1).name = 'par1';
Parameters(1).values = {'a','b','c'};

Parameters(2).name = 'par2';
Parameters(2).values = {'d','e'};

output = prepvalidation(Parameters);

err = ~isequal(size(output),[2*3 2]);
data = [];
maxerr = [];

end