function [pass,maxerr] = test(opt)

% Check that fitparamodel() can fit signals with non-square kernel matrices

t = linspace(0,2,100);
r = linspace(1,6,200);
parain = [3,0.5];
P = rd_onegaussian(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;
par0 = [2 0.1];
[parafit,Pfit] = fitparamodel(S,@rd_onegaussian,r,K,par0);

% Pass 1: distirbution is well fitted
pass(1) = all(abs(Pfit - P) < 1e-5);
% Pass 2: fit parameters agree
pass(2) = all(abs(parafit - parain) < 1e-3);
% Pass 3: the dimensions are correct
pass(3)  = length(Pfit) > length(S);

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