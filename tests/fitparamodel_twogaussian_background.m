function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a two-Gaussian model with background

t = linspace(0,5,300);
r = linspace(2,6,300);
parIn = [2.5 0.5 4 0.5 0.4];
P = dd_twogauss(r,parIn);
K = dipolarkernel(t,r);
B = bg_exp(t,0.15);
lam = 0.25;
V = (1 - lam + lam*K*P).*B;
KB = dipolarkernel(t,r,lam,B);
par0 = [2 0.1 5 0.1 0.1];

[~,Pfit] = fitparamodel(V,@dd_twogauss,r,KB,par0);

%Pass: distance distribution is well fitted
pass = all(abs(Pfit - P) < 1e-5);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end