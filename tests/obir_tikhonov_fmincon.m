function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
Distribution = rd_twogaussian(r,[2,0.3,3.5,0.3,0.5]);

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
RegParam = 40;
Result = obir(S,K,r,'tikhonov',RegMatrix,RegParam,'NoiseLevelAim',NoiseLevel,'Solver','fmincon','axishandle',axhandle);

RegResult = fitregmodel(S,K,r,RegMatrix,'tikhonov',RegParam);

err = norm(Result - Distribution) > norm(RegResult - Distribution);
maxerr = norm(Result - Distribution);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,Distribution,'k') 
    plot(r,Result,'b')
    plot(r,RegResult,'r')
    legend('truth','OBIR','Tikh')
end

end