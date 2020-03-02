function [pass,maxerr] = test(data,opts)

% Check that prepvalidation's random permutation works 

currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private

Parameters.par1 = [1 2 3];
Parameters.par2 = [1 2 3];
out1 = prepvalidation(Parameters,'randperm',true);
out2 = prepvalidation(Parameters,'randperm',false);

% Pass 1: random permutation works
pass(1) = ~isequal(out1,out2);
% Pass 2: the output dimensions are correct
pass(2) = isequal(size(out1),[3*3 2]);
% Pass 3: the output inner dimensions are correct
pass(3) = ~isequal([out1{1,1} out1{2,1} out1{3,1}],[1 2 3]);

pass = all(pass);
 
maxerr = NaN;

cd(currentpath)

end