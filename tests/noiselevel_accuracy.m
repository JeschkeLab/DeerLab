function [err,data,maxerr] = test(opt,olddata)

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.5];
Distribution = rd_onegaussian(r,InputParam);
Distribution = Distribution/sum(Distribution);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

rng(2)
Noise = rand(Dimension,1);
Noise = Noise - mean(Noise);
Noise = 0.02*Noise/Noise(1);

DipEvoFcn = DipEvoFcn + Noise;

truelevel = std(Noise);
approxlevel = noiselevel(DipEvoFcn);

err = abs(approxlevel - truelevel)>1e-2;
maxerr = abs(approxlevel - truelevel);
data = [];


end
