function [err,data,maxerr] = test(opt,olddata)

S1 = ones(40,1);
S2 = ones(70,1);
S3 = ones(120,1);
S = {S1,S2,S3};
NoiseLevel = [0.1 0.05 0.09];

N = cellfun(@length,S);
wref = 1./(N.*NoiseLevel);
wref = wref/sum(wref);

w = globalweights(S,NoiseLevel);

maxerr = max(abs(w-wref));
err = maxerr>1e-6;
data = [];

end
