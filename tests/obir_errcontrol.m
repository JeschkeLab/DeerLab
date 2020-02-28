function [pass,maxerr] = test(opt)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.3]);
rng(2)
K = dipolarkernel(t,r);
S = dipolarsignal(t,r,P,'noiselevel',0.05,'phase',pi/2);

alpha = 0.5;

try
    Pfit = obir(S,K,r,'tikh',alpha);
    err = true;
catch
   err = false;
end

maxerr = 0;
 

end