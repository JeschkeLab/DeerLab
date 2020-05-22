function [pass,maxerr] = test(opt)

% Test Tikhonov regularization with the fnnls solver (free)

rng(1)
t = linspace(0,3,200);
r = linspace(2,4,200);
P = dd_gauss(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.02);
[Pfit,Pci] = fitregmodel(S,K,r,'tikhonov','aic');

%Pass : fnnls manages to fit the distribution
pass = all(abs(Pfit - P) < 3e-1);

maxerr = max(abs(Pfit - P));

if opt.Display
   hold on
   fill([r fliplr(r)],[Pci(1,:) flipud(Pci(2,:))],'r','linestyle','none','facealpha',0.4)
   plot(r,P,'k',r,Pfit,'r')
   hold off
   legend('95% CI','truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end