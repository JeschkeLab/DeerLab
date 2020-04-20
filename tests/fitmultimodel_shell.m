function [pass,maxerr] = test(opt)

% Test fitmultimodel() using a model without built-in multi-component models

t = linspace(0,2,100);
r = linspace(0,7,300);
P = dd_shell(r,[1.5 0.5]) + dd_shell(r,[2.5 0.6]);
P = P/sum(P)/mean(diff(r));
K = dipolarkernel(t,r);
S = K*P;
[Pfit,~,~,~,Nopt,~,~] = fitmultimodel(S,K,r,@dd_shell,4,'aicc','TolFun',1e-5);

% Pass 1: distribution is well fitted
pass(1) = all(abs(Pfit - P) < 7e-2);
% Pass 2: optimal number of components identified correctly
pass(2) = Nopt==2;

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