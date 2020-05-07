function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a three Gaussian model

t = linspace(0,5,300);
r = linspace(2,6,300);
parIn = [3 0.3 0.3 4 0.3 0.3 5 0.3];
P = dd_gauss3(r,parIn);
K = dipolarkernel(t,r);
V = K*P;

rng(67)
lb = parIn - [0.5 0.2 0.2 1 0.2 0.2 0.5 0.2];
ub = parIn + [0.5 0.2 0.2 1 0.2 0.2 1   0.2];
[~,FitP] = fitparamodel(V,@dd_gauss3,r,K,'multistart',5,'Lower',lb,'Upper',ub);

% Pass: distance distribution is well fitted
pass = all(abs(FitP - P) < 1e-10);

maxerr = max(abs(FitP - P));
 
if opt.Display
   plot(r,P,'k',r,FitP,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end