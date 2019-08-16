function [err,data,maxerr] = test(opt,olddata)

currentpath = pwd;
cd(fileparts(mfilename('fullpath')))
cd ../functions/private


cachedData =  java.util.LinkedHashMap;

for i=1:25
hashKey = datahash(i);
Output = {i};
cachedData = addcache(cachedData,hashKey,Output);
end

keys = 6:25;
cacheshould = keys;
for i = 1:cachedData.size
    hashKey = datahash(keys(i));
    out = java2mat({cachedData.get(hashKey)});
    cacheis(i) = out{:};
end

err(1) = cachedData.size > 20;
err(2) = any(cacheis ~= cacheshould);
err = any(err);
data = [];
maxerr = [];

cd(currentpath)

end