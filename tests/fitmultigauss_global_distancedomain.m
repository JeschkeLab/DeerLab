function [pass,maxerr] = test(opt)

% Test that fitmultigauss works with multiple signals using global fit

r = linspace(2,6,300);
parIn = [4 0.2 4 1 3 0.4 0.4 0.4];
P = dd_gauss3(r,parIn);

t1 = linspace(0,2,200);
K1 = dipolarkernel(t1,r);
S1 = K1*P;

t2 = linspace(0,5,200);
K2 = dipolarkernel(t2,r);
S2 = K2*P;

t3 = linspace(0,7,200);
K3 = dipolarkernel(t3,r);
S3 = K3*P;

Pfit = fitmultigauss({S1,S2,S3},{K1,K2,K3},r,5,'aicc');

% Pass: distribution is well fitted
pass = all(abs(Pfit - P) < 7e-2);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end