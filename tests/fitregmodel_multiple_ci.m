function [pass,maxerr] = test(opt)

% Test that multiple confidence levels can be requested via options

t = linspace(0,8,300);
r = linspace(1,8,300);
parIn  = [4.5,0.5];
P = dd_gauss(r,parIn);

K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.05);

[Pfit,Pci] = fitregmodel(S,K,r,'tikh','aic','confidencelevel',[0.5 0.95]);

Pci50 = Pci{1};
Pci95 = Pci{2};

% Pass: confidence intervals behave as expected
pass = all(all(abs(Pfit - Pci50) <= abs(Pfit - Pci95)));

maxerr = NaN;

if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   hold on
   fill([r fliplr(r)],[Pci50(:,1); flipud(Pci50(:,2))],'r','FaceAlpha',0.5,'LineStyle','none')
   fill([r fliplr(r)],[Pci95(:,1); flipud(Pci95(:,2))],'r','FaceAlpha',0.25,'LineStyle','none')
   hold off
   legend('truth','fit','50% CI','95% CI')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end