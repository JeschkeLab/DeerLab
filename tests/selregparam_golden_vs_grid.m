function [pass,maxerr] = test(opt)

N = 100;
dt = 0.008;
t = linspace(0,dt*N,N);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
S = dipolarsignal(t,r,P,'noiselevel',0.05);

tic
[alpha1] = selregparam(S,K,r,'tikhonov','aic','Search','golden');
tictoc1 = toc;
tic
[alpha2] = selregparam(S,K,r,'tikhonov','aic','Search','grid');
tictoc2 = toc;

err(1) = (alpha1 - alpha2) > 1e-1;
err(2) = tictoc1 > tictoc2;

pass = all(err);
 
maxerr = max(abs(alpha1 - alpha2));


end