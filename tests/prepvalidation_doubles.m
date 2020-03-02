function [pass,maxerr] = test(data,opts)

% Check that prepvalidation() works with doubles

currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private

Parameters.par1 = linspace(1,50,10);
Parameters.par2 = linspace(50,100,10);
output = prepvalidation(Parameters);

% Pass: the output dimensions are correct
pass = isequal(size(output),[10*10 2]);
 
maxerr = NaN;

cd(currentpath)


end