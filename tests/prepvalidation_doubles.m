function [err,data,maxerr] = test(data,opts)

Parameters.par1 = linspace(1,50,10);

Parameters.par2 = linspace(50,100,10);

output = prepvalidation(Parameters);

err = ~isequal(size(output),[10*10 2]);
data = [];
maxerr = [];

end