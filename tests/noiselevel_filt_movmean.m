function [pass,maxerr] = test(opt)

% Test noiselevel() of complex signals

rng(1)
t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_gauss(r,[3 0.5]);
lam = 0.25;
B = bg_exp(t,1.5);
V = dipolarsignal(t,r,P,lam,B,'noiselevel',0.03);

noise = whitegaussnoise(t,0.03);

truelevel = std(noise);
approxlevel = noiselevel(V,'movmean');

% Pass: the noise level is well estimated
pass = abs(approxlevel - truelevel) < 1e-2;

maxerr = abs(approxlevel - truelevel);
 
end
