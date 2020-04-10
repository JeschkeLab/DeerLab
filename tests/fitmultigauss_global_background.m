function [pass,maxerr] = test(opt)

% Check that fitmultigauss do global fit with different backgrounds

rng(1)

r = linspace(2,6,200);
parIn = [4 0.2, 3.5 1 0.4];
P = dd_gauss2(r,parIn);

t1 = linspace(0,5,100);
lam1 = 0.2;
k1 = 0.2;
B1 = bg_exp(t1,k1);
V1 = dipolarsignal(t1,r,P,lam1,B1);

t2 = linspace(0,3,100);
lam2 = 0.3;
k2 = 0.3;
B2 = bg_exp(t2,k2);
V2 = dipolarsignal(t2,r,P,lam2,B2);

[Pfit,parfit] = fitmultigauss({V1,V2},{t1,t2},r,3,'aicc','background',@bg_exp);

% Pass 1: distribution is well fitted
pass(1) = all(abs(Pfit - P) < 8e-1);
% Pass 2-3: modulation depth and background parameters of individual signals are well fitted
pass(2) = all(abs(parfit(9:10) - [lam1 k1]) < 1e-2);
pass(3) = all(abs(parfit(11:12) - [lam2 k2]) < 1e-2);

maxerr = max(abs(Pfit - P));

if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end