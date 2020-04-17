function [pass,maxerr] = test(opt)

% Check that fitmultigauss global fitting with boundaries
rng(1)

r = linspace(2,6,200);
parIn = [4 0.2, 3.5 1 0.4];
P = dd_gauss2(r,parIn);

t1 = linspace(0,5,100);
V1 = dipolarsignal(t1,r,P);

t2 = linspace(0,3,100);
V2 = dipolarsignal(t2,r,P);

Plower = [1 0.2];
Pupper = [20 5];

[Pfit] = fitmultigauss({V1,V2},{t1,t2},r,3,'aicc','upper',Pupper,'lower',Plower,'tolfun',1e-4);

% Pass 1: distribution is well fitted
pass = all(abs(Pfit - P) < 8e-1);

maxerr = max(abs(Pfit - P));

if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end