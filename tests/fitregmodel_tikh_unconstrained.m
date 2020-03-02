function [pass,maxerr] = test(opt)

% Test unconstrained Tikhonov regularization

t = linspace(0,2,200);
r = linspace(2,5,100);
P = rd_onegaussian(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;
alpha = 1;

Pfit = fitregmodel(S,K,r,'tikhonov',alpha,'NonNegConstrained',false);

pass = all(abs(Pfit - P) < 4e-1);

maxerr = max(abs(Pfit - P));

if opt.Display
   plot(r,P,'k',r,Pfit)
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end