function [pass,maxerr] = test(opt)

% Check error control of nosielevel() towards wrong inputs

rng(1)
t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_onegaussian(r,[3 0.5]);
K = dipolarkernel(t,r);
S = K*P;

noise = rand(numel(t),1);
noise = noise - mean(noise);
noise = 0.02*noise/noise(1);
S = S + noise;

% Pass: asking to remove too many points
try
    noiselevel(S,1000);
    pass = false;
catch
    pass = true; 
end

maxerr = NaN;
 


end
