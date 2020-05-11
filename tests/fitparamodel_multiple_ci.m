function [pass,maxerr] = test(opt)

% Test that multiple confidence levels can be requested via options

t = linspace(0,8,300);
r = linspace(1,8,300);
parIn  = [4.5,0.5];
P = dd_gauss(r,parIn);

K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.05);

par0 = [2 0.2];
[parfit,Pfit,parCI,Pci] = fitparamodel(S,@dd_gauss,r,K,par0,'confidencelevel',[0.5 0.95]);

parci1 = parCI{1};
parci2 = parCI{2};

% Pass 1-2: confidence intervals behave as expected
pass = all(all(abs(parfit - parci1.') < abs(parfit - parci2.')));

maxerr = NaN;

if opt.Display
   plot(r,P,'k',r,Pfit,'r',r,Pci{2},'r--')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end