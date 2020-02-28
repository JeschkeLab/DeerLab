function [pass,maxerr] = test(opt)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_twogaussian(r,[2,0.3,3.5,0.3,0.5]);

K = dipolarkernel(t,r);
RegMatrix =  regoperator(Dimension,2);
DipEvoFcn = K*P;
NoiseLevel = 0.05;
Noise = whitegaussnoise(Dimension,NoiseLevel);
S = DipEvoFcn+Noise;

%Set optimal regularization parameter (found numerically lambda=0.13)
OptParam = 12;
OptHuber = 1.35;

if opt.Display
    figure(8),clf
    axhandle = axes();
else
    axhandle = [];
end

Result = obir(S,K,r,'huber',OptParam,'DivergenceStop',true,'NoiseLevelAim',NoiseLevel,'Solver','fnnls','Huberparam',OptHuber,'axishandle',axhandle);

RegResult = fitregmodel(S,K,r,'huber',OptParam);

err = norm(Result - P) > norm(RegResult - P);
maxerr = norm(Result - P);
 

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Result,'b')
    plot(r,RegResult,'r')
    legend('truth','OBIR','Huber')
end

end