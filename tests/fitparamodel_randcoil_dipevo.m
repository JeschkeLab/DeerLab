function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a random coil model

t = linspace(0,3,400);
r = linspace(1,6,150);
parIn = [20,0.2,0.49];
P = dd_randcoil(r,parIn);
K = dipolarkernel(t,r);
S = K*P;

[parFit,Pfit] = fitparamodel(S,@dd_randcoil,r,K);

%Pass 1: distance distribution is well fitted
pass(1) = all(abs(Pfit - P) < 1e-7);

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
