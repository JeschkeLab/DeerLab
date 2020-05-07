function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a wormchain model
t = linspace(0,8,300);
r = linspace(1,8,300);
parIn  = [4.5,0.5];
P = dd_gauss(r,parIn);

K = dipolarkernel(t,r);
S = K*P;

par0 = [2 0.2];
[parFit,Pfit,parCI] = fitparamodel(S,@dd_gauss,r,K,par0,'solver','lsqnonlin');

% Pass 1: distance distribution is well fitted
pass(1) = all(abs(Pfit - P) < 1e-5);
% Pass 2: model parameters are well fitted
pass(2) = all(abs(parFit - parIn) < 1e-3);
% Pass 3: confidence intervals are exactly zero
pass(3) = all(all(abs(parFit - parCI.') < 1e-10));

maxerr = max(abs(Pfit - P));

 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end