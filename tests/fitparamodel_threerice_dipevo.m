function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a three Rician model

t = linspace(0,5,300);
r = linspace(2,6,300);
parIn = [2 0.4 0.3 3.5 0.3 0.3 5 0.3];
P = dd_rice3(r,parIn);
K = dipolarkernel(t,r);
S = K*P;

[~,FitP] = fitparamodel(S,@dd_rice3,r,K,0.75*parIn);

%Pass: distance distribution is well fitted
pass = all(abs(FitP - P) < 1e-2);

maxerr = max(abs(FitP - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end