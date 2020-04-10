function [pass,maxerr] = test(opt)

% Check that fitmultigauss global fitting with background boundaries
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

Plower = [1 0.2];
Pupper = [20 5];
B1lower = [0.1 0.1];
B1upper = [0.4 0.4];
B2lower = [0.1 0.1];
B2upper = [0.4 0.4];
upper = [Pupper B1upper B2upper];
lower = [Plower B1lower B2lower];

[Pfit,parfit] = fitmultigauss({V1,V2},{t1,t2},r,3,'aicc','background',@bg_exp,'upper',upper,'lower',lower,'tolfun',1e-4);

% % Pass 1: distribution is well fitted
pass(1) = all(abs(Pfit - P) < 8e-1);
% Pass 2-3: modulation depth and background parameters of individual signals are well fitted
pass(2) = all(abs(parfit(end-3:end-2) - [lam1 k1]) < 1e-2);
pass(3) = all(abs(parfit(end-1:end) - [lam2 k2]) < 1e-2);

maxerr = max(abs(Pfit - P));

if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end