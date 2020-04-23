function [pass,maxerr] = test(opt)

% Check that fitmultimodel() works with time-axis as input

t = linspace(0,5,100);
r = linspace(2,6,300);
InputParam = [4 0.2 0.4 4 1 0.4 3 0.4];
P = dd_gauss3(r,InputParam);
V = dipolarsignal(t,r,P);
[Pfit,~,Pci] = fitmultimodel(V,t,r,@dd_gauss,5,'aicc');

% Pass: distribution is well fitted
pass = all(abs(Pfit - P) < 8e-1);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   hold on
   plot(r,P,'k',r,Pfit,'r')
   fill([r fliplr(r)], [Pci(:,1); flipud(Pci(:,2))],'r','Linestyle','none','facealpha',0.25)
   hold off
   legend('truth','fit','ci')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end