function [pass,maxerr] = test(opt)

% Check that fitmultimodel() works well with upper/lower bounds (3 parameters)

t = linspace(0,3,100);
r = linspace(0,6,300);
P = dd_gengauss(r,[2.5 0.5 5]) + 0.8*dd_gengauss(r,[3 0.7 2]);
P = P/sum(P)/mean(diff(r));
K = dipolarkernel(t,r);
S = K*P;

Pfit = fitmultimodel(S,K,r,@dd_gengauss,3,'aic','Lower',[1 0.2 0.25],'Upper',[20 5 15]);

% Pass 1-3: all distributions are well fitted
pass = any(abs(Pfit - P) < 7e-5);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit)
   legend('truth','fit 1','fit 2','fit 3')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end