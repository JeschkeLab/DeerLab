function [pass,maxerr] = test(opt)

% Check that fitmultigauss do global fit with different backgrounds

rng(1)

r = linspace(2,6,200);
parIn = [4 0.2 4 1 3 0.4 0.4 0.4];
P = dd_gauss3(r,parIn);

t1 = linspace(0,5,100);
B1 = bg_exp(t1,0.2);
V1 = dipolarsignal(t1,r,P,'moddepth',0.2,'background',B1);

t2 = linspace(0,3,100);
B2 = bg_exp(t2,0.3);
V2 = dipolarsignal(t2,r,P,'moddepth',0.3,'background',B2);

[Pfit,parfit] = fitmultigauss({V1,V2},{t1,t2},r,5,'aicc','background',@bg_exp);

% Pass 1: distribution is well fitted
pass(1) = all(abs(Pfit - P) < 8e-1);
% Pass 2-3: modulation depth and background parameters of individual signals are well fitted
pass(2) = all(abs(parfit(6:7) - 0.2) < 5e-2);
pass(3) = all(abs(parfit(8:9) - 0.3) < 5e-2);

maxerr = max(abs(Pfit - P));

if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end