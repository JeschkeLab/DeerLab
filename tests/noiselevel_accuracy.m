function [pass,maxerr] = test(opt)

% Check that noiselevel() estimates the noise level accurately

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

truelevel = std(noise);
approxlevel = noiselevel(S);

% Pass: the noise level is well estimated
pass = abs(approxlevel - truelevel) < 1e-2;

maxerr = abs(approxlevel - truelevel);
 


end
