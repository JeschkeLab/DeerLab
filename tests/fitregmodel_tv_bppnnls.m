function [pass,maxerr] = test(opt)

% Test TV regularization with the bppnnls solver (free)

rng(1)
t = linspace(0,3,200);
r = linspace(1,5,100);
P = dd_onegauss(r,[3,0.2]);
K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.01);
alpha = 0.01356;

Pfit = fitregmodel(S,K,r,'tv',alpha,'Solver','bppnnls');

error = abs(Pfit - P);

%Pass : bppnnls manages to fit the distribution
pass = all(error < 3e-1);

maxerr = max(abs(Pfit - P));

if opt.Display
   plot(r,P,'k',r,Pfit)
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end