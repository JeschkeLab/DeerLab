function [pass,maxerr] = test(opt)

% Test a distance-domain fit of a random coil model

t = linspace(0,3,400);
r = linspace(1,6,150);
N = 20; % number of residues
nu = 0.49; % scaling exponent
parIn = [N,nu];
P = rd_randcoil(r,parIn);
K = dipolarkernel(t,r);
S = K*P;

[parFit,Pfit] = fitparamodel(S,@rd_randcoil,r,K);

%Pass 1: distance distribution is well fitted
pass(1) = all(abs(Pfit - P) < 1e-7);
%Pass 2: model parameters are well fitted
pass(2) = all(abs(parFit - parIn)./parIn < 5e-1);

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
