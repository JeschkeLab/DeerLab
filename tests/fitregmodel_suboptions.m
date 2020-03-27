function [pass,maxerr] = test(opt)

% Test that selregparam options can be passed via fitregmodel

rng(1)
t = linspace(0,3,200);
r = linspace(1,5,100);
P = dd_onegauss(r,[3,0.2]);
K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.01);
alpha = 0.01356;

Pfit = fitregmodel(S,K,r,'tv',alpha,'Solver','bppnnls','Search','grid');

error = abs(Pfit - P);

%Pass : fitregmodel manages to pass the options to selregparam
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