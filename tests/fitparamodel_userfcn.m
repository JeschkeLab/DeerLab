function [pass,maxerr] = test(opt)

% Test a distance-domain fit of custom defined model

t = linspace(0,5,200);
r = linspace(1,6,300);
P = dd_onegauss(r,[3,0.2]);
K = dipolarkernel(t,r);
S = K*P;

InitialGuess = [3.5 0.3];
fcnhandle = @(r,p)exp(-((r - p(1))/(p(2))).^2);

[~,Pfit] = fitparamodel(S,fcnhandle,r,K,InitialGuess);

%Pass 1: distance distribution is well fitted
pass(1) = all(abs(Pfit - P) < 1e-9);
%Pass 2: dimensions are right
pass(2)  = length(Pfit) > length(S);

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