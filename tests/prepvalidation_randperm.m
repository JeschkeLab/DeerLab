function [err,data,maxerr] = test(data,opts)

currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private

Parameters.par1 = [1 2 3];
Parameters.par2 = [1 2 3];

output = prepvalidation(Parameters,'randperm',true);

err(1) = ~isequal(size(output),[3*3 2]);
err(2) = isequal([output{1,1} output{2,1} output{3,1}],[1 2 3]);
err = any(err);
data = [];
maxerr = NaN;

cd(currentpath)


end