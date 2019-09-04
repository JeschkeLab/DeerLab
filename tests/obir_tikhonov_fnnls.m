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
P = 0.5*P1 + 0.5*P2;
P = P/sum(P)/mean(diff(r));

K = dipolarkernel(t,r);
RegMatrix =  regoperator(Dimension,2);
DipEvoFcn = K*P;
NoiseLevel = 0.05;
Noise = whitegaussnoise(Dimension,NoiseLevel);
S = DipEvoFcn+Noise;

%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 40;

if opt.Display
    figure(8),clf
    axhandle = plot(r,NaN*P);
else
    axhandle = [];
end

Result = obir(S,K,r,'tikhonov',RegMatrix,RegParam,'NoiseLevelAim',NoiseLevel,'Solver','fnnls','Axishandle',axhandle);
RegResult = fitregmodel(S,K,r,RegMatrix,'tikhonov',RegParam);

err = norm(Result - P) > norm(RegResult - P);
maxerr = norm(Result - P);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Result,'b')
    plot(r,RegResult,'r')
    legend('truth','OBIR','Tikh')
end

end