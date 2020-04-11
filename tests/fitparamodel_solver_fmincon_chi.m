function [pass,maxerr] = test(opt)

% Check that fitparamodel works with fmincon (toolbox) solver using the Chi^2 cost function 

rng(2)
t = linspace(0,3,200);
r = linspace(2,6,150);
parIn = [3,0.3];
P = dd_gauss(r,parIn);
K = dipolarkernel(t,r);
S = dipolarsignal(t,r,P,'noiselevel',0.01);

[parFit,Pfit] = fitparamodel(S,@dd_gauss,r,K,'costmodel','chi2red');

%Pass 1-2: fmincon finds the correct solution with the Chi^2
pass(1) = any(abs(Pfit - P) < 1e-1);
pass(2) = all(abs(parFit - parIn) < 1e-1);

pass = all(pass);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end