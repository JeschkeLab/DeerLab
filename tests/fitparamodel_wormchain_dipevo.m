function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a wormchain model

t = linspace(0,5,300);
r = linspace(1,6,300);
parIn = [4 4.5];
P = dd_wormchain(r,[2.7,20]);

K = dipolarkernel(t,r);
S = K*P;

[~,Pfit] = fitparamodel(S,@dd_wormchain,r,K,'Solver','lsqnonlin');

%Pass 1: distance distribution is well fitted
pass = any(abs(Pfit - P) < 5e-3);

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