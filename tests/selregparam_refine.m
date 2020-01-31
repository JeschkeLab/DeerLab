function [err,data,maxerr] = test(opt,olddata)
Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);
P = P/sum(P);

K = dipolarkernel(t,r);
RegMatrix = regoperator(Dimension,2);
DipEvoFcn = K*P;


Noise = rand(Dimension,1);
Noise = Noise - mean(Noise);
Noise = 0.05*Noise/max(Noise);

RegParamSet = regparamrange(K,RegMatrix);
goal = 0.01;
[~,Functionals,RegParams] = selregparam(DipEvoFcn + Noise,K,r,'tikhonov',{'aic','gcv'},'Refine',true);

err = length(RegParams) == length(RegParamSet);
data = [];
maxerr = NaN;

if opt.Display
   figure(8),clf
hold on
a = Functionals{1};
plot((RegParams(1:end-60)),a(1:end-60),'.')
plot((RegParams(end-60:end)),a(end-60:end),'.')
end

end