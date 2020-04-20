function [pass,maxerr] = test(opt)

% Test that fitmultimodel() works with multiple signals using global fit

r = linspace(2,6,300);
parIn = [4 0.2 4 1 3 0.4 0.4 0.4];
P = dd_gauss3(r,parIn);

t1 = linspace(0,2,200);
K1 = dipolarkernel(t1,r);
S1 = K1*P;

t2 = linspace(0,5,200);
K2 = dipolarkernel(t2,r);
S2 = K2*P;

Pfit = fitmultimodel({S1,S2},{K1,K2},r,@dd_gauss,4,'aicc');

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