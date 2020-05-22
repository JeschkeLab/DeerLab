function [pass,maxerr] = test(opt)

% Test Huber regularization using the fmincon solver (toolbox)

t = linspace(0,4,300);
r = linspace(2,6,100);
P = dd_gauss(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;

Pfit = fitregmodel(S,K,r,'huber','aic','Solver','fmincon','HuberParam',1.35);


error = abs(Pfit - P);
% Pass: the distribution is well fitted
pass = all(error < 2e-1);

maxerr = max(error);
 

if opt.Display
   plot(r,P,r,Pfit)
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end