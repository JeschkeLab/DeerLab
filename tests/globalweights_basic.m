function [pass,maxerr] = test(opt)

% Check basic functionality of the globalweights() function

S1 = ones(40,1);
S2 = ones(70,1);
S3 = ones(120,1);
S = {S1,S2,S3};
levels = [0.1 0.05 0.09];

N = cellfun(@length,S);
wref = 1./(N.*levels);
wref = wref.'/sum(wref);

w = globalweights(S,levels);

maxerr = max(abs(w - wref));

% Pass: the function output matches the refernce
pass = maxerr < 1e-6;
 

end
