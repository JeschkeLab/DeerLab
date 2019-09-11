function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
NoiseLevel = 0.01;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);
rng(2)
K = dipolarkernel(t,r);
DipEvoFcn = K*P;
Noise = whitegaussnoise(Dimension,NoiseLevel);
S = DipEvoFcn + Noise;
rng('default')
RegParam = 0.6310;
TikhResult1 = fitregmodel(S,K,r,'huber',RegParam,'Solver','fnnls');
err(1) = any(abs(TikhResult1 - P)>9e-2);
maxerr = max(abs(TikhResult1 - P));

err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    subplot(121)
    hold on
    plot(t,S)
    plot(t,K*TikhResult1)
    subplot(122)
    hold on
    plot(r,P,'k') 
    plot(r,TikhResult1,'r')
end

end