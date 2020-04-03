function [pass,maxerr] = test(opt)

% Check that noiselevel() accepts a manual point for the estimation

rng(1)
t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_gauss(r,[3 0.5]);
K = dipolarkernel(t,r);
S = K*P;

noise = rand(numel(t),1);
noise = noise - mean(noise);
noise = 0.02*noise/noise(1);
S = S + noise;

approxlevel1 = noiselevel(S);
approxlevel2 = noiselevel(S,1/5*numel(t));

% Pass: the noise level is well estimated
pass = abs(approxlevel1 - approxlevel2) < 1e-10;

maxerr = abs(approxlevel1 - approxlevel2);
 


end
