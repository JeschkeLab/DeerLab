function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;
Noise = whitegaussnoise(Dimension,0.02);
%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.1;
RegMatrix = regoperator(Dimension,3);
Resultfnnls = fitregmodel(DipEvoFcn+Noise,K,r,RegMatrix,'tv',RegParam,'Solver','fnnls');

err = any(abs(Resultfnnls - P)>9e-2);

maxerr = max(abs(Resultfnnls - P));

data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k')
    plot(r,Resultfnnls,'b') 
    axis tight
end

end