function [pass,maxerr] = test(opt)

% Test fitmultimodel() using a model without built-in multi-component
% models and more than two parameters

t = linspace(0,2,100);
r = linspace(0,7,300);
P = dd_gengauss(r,[2.5 0.5 5]) + 0.8*dd_gengauss(r,[3 0.7 2]);
P = P/sum(P)/mean(diff(r));
K = dipolarkernel(t,r);
S = K*P;
Pfit = fitmultimodel(S,K,r,@dd_gengauss,4,'aicc','TolFun',1e-5);

% Pass: distribution is well fitted
pass = all(abs(Pfit - P) < 4e-1);

maxerr = max(abs(Pfit - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end