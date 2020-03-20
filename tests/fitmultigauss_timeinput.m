function [pass,maxerr] = test(opt)

% Check that fitmultigauss() works with time-axis as input

t = linspace(0,5,100);
r = time2dist(t);
InputParam = [4 0.2 4 1 3 0.4 0.4 0.4];
P = dd_threegauss(r,InputParam);
V = dipolarsignal(t,r,P);
Pfit = fitmultigauss(V,t,r,5,'aicc');

% Pass: distribution is well fitted
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