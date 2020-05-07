function [pass,maxerr] = test(opt)

% Test Tikhonov regularization with the fnnls solver (free)

rng(1)
r = linspace(2,4,200);
P = dd_gauss(r,[3,0.5]);

t1 = linspace(0,3,200);
K1 = dipolarkernel(t1,r);
S1 = K1*P + whitegaussnoise(t1,0.02);

t2 = linspace(0,6,200);
K2 = dipolarkernel(t2,r);
S2 = K2*P + whitegaussnoise(t2,0.03);

alpha = 2;

[Pfit,Pci] = fitregmodel({S1,S2},{K1,K2},r,'tikhonov',alpha,'confidencelevel',0.90);

%Pass : fnnls manages to fit the distribution
pass = all(abs(Pfit - P) < 3e-1);

maxerr = max(abs(Pfit - P));

if opt.Display
   hold on
   fill([r fliplr(r)],[Pci(1,:) flipud(Pci(2,:))],'r','linestyle','none','facealpha',0.4)
   plot(r,P,'k',r,Pfit,'r')
   hold off
   legend('90% CI','truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end