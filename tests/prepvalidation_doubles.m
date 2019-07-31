function [err,data,maxerr] = test(data,opts)

Parameters(1).name = 'par1';
Parameters(1).values = linspace(1,50,10);

Parameters(2).name = 'par2';
Parameters(2).values = linspace(50,100,10);

output = prepvalidation(Parameters);

err = ~isequal(size(output),[10*10 2]);
data = [];
maxerr = [];

end