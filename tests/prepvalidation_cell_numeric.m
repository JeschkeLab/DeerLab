function [pass,maxerr] = test(data,opts)

currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private


Parameters.par1 = num2cell(rand(5,100),2);
Parameters.par2 = num2cell(rand(6,100),2);

output = prepvalidation(Parameters);

err(1) = ~isequal(size(output),[6*5 2]);
err(2) = length(output{1,1})~= 100;
pass = all(err);
 
maxerr = NaN;

cd(currentpath)


end