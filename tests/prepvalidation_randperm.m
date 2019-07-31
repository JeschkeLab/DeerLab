function [err,data,maxerr] = test(data,opts)

Parameters(1).name = 'par1';
Parameters(1).values = [1 2 3];

Parameters(2).name = 'par2';
Parameters(2).values = [1 2 3];

output = prepvalidation(Parameters,'randperm',true);

err(1) = ~isequal(size(output),[3*3 2]);
err(2) = isequal([output{1,1} output{2,1} output{3,1}],[1 2 3]);
err = any(err);
data = [];
maxerr = [];

end