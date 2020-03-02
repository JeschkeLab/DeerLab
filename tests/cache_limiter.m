function [pass,maxerr] = test(opt)

% Check that the Java cache memory is properly limited 

%Temporarely move to private directory
currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private

%Create Java linked hash map
cachedData =  java.util.LinkedHashMap;
%Overfill hash map
for i=1:25
hashKey = datahash(i);
Output = {i};
cachedData = addcache(cachedData,hashKey,Output);
end
%Retrieve chached data
keys = 6:25;
cacheshould = keys;
for i = 1:cachedData.size
    hashKey = datahash(keys(i));
    out = java2mat({cachedData.get(hashKey)});
    cacheis(i) = out{:};
end

% Pass 1: cache size is correctly limited
pass(1) = cachedData.size <= 20;
% Pass 2: all cached data is the one it should be
pass(2) = all(cacheis == cacheshould);

pass = all(pass);
 
maxerr = NaN;

%Return to original directory
cd(currentpath)

end