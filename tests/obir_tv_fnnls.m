function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
dt = 0.008;
rng(2)
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P1 = rd_onegaussian(r,[2,0.3]);
P2 = rd_onegaussian(r,[3.5,0.3]);
P = 0.5*P1 + 0.5*P2;
P = P/sum(P)/mean(diff(r));

K = dipolarkernel(t,r);
DipEvoFcn = K*P;
NoiseLevel = 0.02;
Noise = whitegaussnoise(Dimension,NoiseLevel);
S = DipEvoFcn+Noise;

if opt.Display
    figure(8),clf
    axhandle = axes();
else
    axhandle = [];
end

%Set optimal regularization parameter (found numerically lambda=0.13)
OptParam = 0.002;
Result = obir(S,K,r,'tv',OptParam,'DivergenceStop',true,'NoiseLevelAim',NoiseLevel,'Solver','fnnls','axishandle',axhandle);
RegResult = fitregmodel(S,K,r,'tv',OptParam);

err = norm(Result - P) > norm(RegResult - P);
maxerr = norm(Result - P);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Result,'b')
    plot(r,RegResult,'r')
    legend('truth','OBIR','TV')
end

end