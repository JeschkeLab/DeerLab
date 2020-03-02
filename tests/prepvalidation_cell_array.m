function [pass,maxerr] = test(data,opts)

% Check that prepvalidation() works with cell arrays containing strings

currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private

Parameters.par1 = {'a','b','c'};
Parameters.par2 = {'d','e'};
output = prepvalidation(Parameters);

% Pass: the output dimensions are correct
pass = isequal(size(output),[6 2]);
 
maxerr = NaN;

cd(currentpath)

end