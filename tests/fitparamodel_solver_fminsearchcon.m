function [pass,maxerr] = test(opt)

% Check that fitparamodel works with fminsearchcon (free) solver 

rng(1)
t = linspace(0,3,200);
r = linspace(2,6,150);
parIn  = [3,0.5];
P = dd_onegauss(r,parIn);

K = dipolarkernel(t,r);
S = K*P;

param0 = [2 0.1];
[parFit,Pfit] = fitparamodel(S,@dd_onegauss,r,K,param0,'Solver','fminsearchcon');

%Pass 1-2: fminsearchcon finds the correct solution
pass(1) = all(abs(Pfit - P) < 3e-5);
pass(2) = all(abs(parFit - parIn) < 1e-3);

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