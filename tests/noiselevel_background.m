function [pass,maxerr] = test(opt)

% Check that noiselevel() estimates the noise level accurately even in the
% presence of a background function

rng(1)
t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_onegaussian(r,[3 0.5]);
K = dipolarkernel(t,r);
lam = 0.25;
B = bg_exp(t,1.5);
V = (1 - lam + lam*K*P).*B;

noise = rand(numel(t),1);
noise = noise - mean(noise);
noise = 0.02*noise/noise(1);
V = V + noise;

truelevel = std(noise);
approxlevel = noiselevel(V);

% Pass: the noise level is well estimated
pass = abs(approxlevel - truelevel) < 1e-2;

maxerr = abs(approxlevel - truelevel);
 
end
