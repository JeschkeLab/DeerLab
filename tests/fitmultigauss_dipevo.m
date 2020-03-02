function [pass,maxerr] = test(opt)

% Test the fit quality of fitmultigauss() with a dipolar evolution function

t = linspace(0,2,100);
r = time2dist(t);
paramtrue = [4 0.2 4 1 3 0.4 0.4 0.4];
P = rd_threegaussian(r,paramtrue);
K = dipolarkernel(t,r);
S = K*P;
Pfit = fitmultigauss(S,K,r,5,'aicc','TolFun',1e-15);

% Pass: distribution is well fitted
pass = all(abs(Pfit - P) < 7e-2);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end