function [pass,maxerr] = test(data,opts)

currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private

Parameters.par1 = {'a','b','c'};
Parameters.par2 = {'d','e'};

output = prepvalidation(Parameters);

err = ~isequal(size(output),[2*3 2]);
 
maxerr = NaN;

cd(currentpath)

end