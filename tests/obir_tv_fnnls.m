function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P1 = rd_onegaussian(r,[2,0.3]);
P2 = rd_onegaussian(r,[3.5,0.3]);
Distribution = 0.5*P1 + 0.5*P2;
Distribution = Distribution/sum(Distribution)/mean(diff(r));

K = dipolarkernel(t,r);
RegMatrix =  regoperator(Dimension,2);
DipEvoFcn = K*Distribution;
NoiseLevel = 0.05;
Noise = whitegaussnoise(Dimension,NoiseLevel);
S = DipEvoFcn+Noise;

if opt.Display
    figure(8),clf
    axhandle = plot(r,NaN*Distribution);
else
    axhandle = [];
end

%Set optimal regularization parameter (found numerically lambda=0.13)
OptParam = 0.1;
Result = obir(S,K,r,'tv',RegMatrix,OptParam,'DivergenceStop',true,'NoiseLevelAim',NoiseLevel,'Solver','fnnls','axishandle',axhandle);
RegResult = fitregmodel(S,K,r,RegMatrix,'tv',OptParam);

err = norm(Result - Distribution) > norm(RegResult - Distribution);
maxerr = norm(Result - Distribution);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,Distribution,'k') 
    plot(r,Result,'b')
    plot(r,RegResult,'r')
    legend('truth','OBIR','TV')
end

end