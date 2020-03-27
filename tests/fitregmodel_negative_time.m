function [pass,maxerr] = test(opt)

% Check that fitregmodel() works with negative times

rng(1)
t = linspace(-2,4,300);
r = linspace(2,6,100);
P = dd_onegauss(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.01);
alpha = 0.31;
Pfit = fitregmodel(S,K,r,'tikhonov',alpha);

error = abs(Pfit - P);
% Pass: the distribution is well fitted
pass = all(error < 5e-2);

maxerr = max(error);
 
if opt.Display
   plot(r,P,r,Pfit)
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end