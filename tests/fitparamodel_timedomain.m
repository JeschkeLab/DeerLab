function [pass,maxerr] = test(opt)

% Test a time-domain fit of a single Gaussian model

t = linspace(0,5,300);
r = linspace(2,6,300);
parIn = [3 0.5];
P = dd_gauss(r,parIn);
K = dipolarkernel(t,r);
V = K*P;
InitialGuess = [2 0.1];

mymodel = @(t,param)K*dd_gauss(r,param);
parFit = fitparamodel(V,mymodel,t,InitialGuess);
Pfit = dd_gauss(r,parFit);


%Pass 1: distance distribution is well fitted
pass(1) = all(abs(Pfit - P) < 1e-5);
%Pass 2: model parameters are well fitted
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