function [pass,maxerr] = test(opt)

% Test Tikhonov regularization with a background-kernel

t = linspace(0,3,200);
r = linspace(2,6,100);
P = dd_onegaussian(r,[3,0.5]);
K = dipolarkernel(t,r);
B = bg_exp(t,0.15);
lam = 0.25;
V = (1 - lam + lam*K*P).*B;
alpha = 0.13;
KB = dipolarkernel(t,r,lam,B);

Pfit = fitregmodel(V,KB,r,'tikhonov',alpha,'Solver','fnnls');

% Pass: the distribution is well fitted
pass = all(abs(Pfit - P) < 3e-1);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   plot(r,P,r,Pfit)
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end