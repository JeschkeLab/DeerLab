function [pass,maxerr] = test(data,opts)

% Check that prepvalidation() works with cell arrays containing numericals

currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private

Parameters.par1 = num2cell(rand(5,100),2);
Parameters.par2 = num2cell(rand(6,100),2);
output = prepvalidation(Parameters);

% Pass 1: the output dimensions are correct
pass(1) = isequal(size(output),[30 2]);
% Pass 1: The dimensions of the array are correct
pass(2) = length(output{1,1}) == 100;

pass = all(pass);
 
maxerr = NaN;

cd(currentpath)


end