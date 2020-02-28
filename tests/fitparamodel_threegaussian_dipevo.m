function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a three Gaussian model

t = linspace(0,5,300);
r = linspace(2,6,300);
parIn = [3 0.3 4 0.3 5 0.3 0.3 0.3];
P = rd_threegaussian(r,parIn);
K = dipolarkernel(t,r);
S = K*P;

[~,FitP] = fitparamodel(S,@rd_threegaussian,r,K);

%Pass: distance distribution is well fitted
pass = all(abs(FitP - P) < 1e-10);

maxerr = max(abs(FitP - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end