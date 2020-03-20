function [pass,maxerr] = test(opt)

%Check that fitparamodel works with negative times

rng(2)
t = linspace(-5,1,300);
r = linspace(2,6,200);
P = dd_onegaussian(r,[4 0.4]);
K = dipolarkernel(t,r);
V = K*P + whitegaussnoise(t,0.02);

[~,Pfit] = fitparamodel(V,@dd_onegaussian,r,K);
error = abs(Pfit - P);

%Pass: solution fits the ground truth
pass = all(error < 7e-2);

maxerr = max(error);

if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end