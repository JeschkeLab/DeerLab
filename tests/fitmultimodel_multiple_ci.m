function [pass,maxerr] = test(opt)

% Test that multiple confidence levels can be requested via options

t = linspace(0,8,300);
r = linspace(1,8,300);
parIn  = [4.5,0.5];
P = dd_gauss(r,parIn);

K = dipolarkernel(t,r);
S = K*P + whitegaussnoise(t,0.05);

[Pfit,parfit,Pci,parci] = fitmultimodel(S,K,r,@dd_gauss,3,'aic','confidencelevel',[0.5 0.95]);

parci50 = parci{1};
parci95 = parci{2};

Pci50 = Pci{1};
Pci95 = Pci{2};

% Pass 1-2: confidence intervals behave as expected
pass(1) = all(all(abs(parfit - parci50.') < abs(parfit - parci95.')));
pass(2) = all(all(abs(Pfit - Pci50) <= abs(Pfit - Pci95)));

pass = all(pass);

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