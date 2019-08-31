function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
NoiseLevel = 0.01;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
Distribution = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;
Noise = whitenoise(Dimension,NoiseLevel);
S = DipEvoFcn + Noise;

RegMatrix = regoperator(Dimension,2);
range = regparamrange(K,RegMatrix);
RegParam = selregparam(range,S,K,RegMatrix,'tikhonov','aicc');

TikhResult1 = fitregmodel(S,K,r,RegMatrix,'tikhonov',RegParam,'Solver','fnnls');
err(1) = any(abs(TikhResult1 - Distribution)>6e-2);
maxerr = max(abs(TikhResult1 - Distribution));

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
    plot(r,Distribution,'k') 
    plot(r,TikhResult1,'r')
end

end