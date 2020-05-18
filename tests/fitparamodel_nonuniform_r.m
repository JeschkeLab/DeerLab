function [pass,maxerr] = test(opt)

%Test distance-domain fit of one Gaussian

t = linspace(0,3,300);
r = sqrt(linspace(1,7^2,200));
parIn  = [3,0.5];
P = dd_gauss(r,parIn);

K = dipolarkernel(t,r);
V = K*P;

par0 = [2 0.1];
[parFit,Pfit] = fitparamodel(V,@dd_gauss,r,K,par0);

%Pass 1: distance distribution is well fitted
pass(1) = all(abs(Pfit - P) < 1e-5);
%Pass 2: model parameters are well fitted
pass(2) = all(abs(parFit - parIn) < 1e-3);

pass = all(pass);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end