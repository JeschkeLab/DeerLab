function [pass,maxerr] = test(opt)

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.5];
P = rd_onegaussian(r,InputParam);
P = P/sum(P);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;
B = exp(-1.5*t)';
V = (DipEvoFcn + 5).*B;
B = B*(1-1/V(1));
V = V/V(1);
V = V./sqrt(B);

rng(2)
Noise = rand(Dimension,1);
Noise = Noise - mean(Noise);
Noise = 0.02*Noise/Noise(1);

V = V + Noise;

truelevel = std(Noise);
approxlevel = noiselevel(V);

err = abs(approxlevel - truelevel)>1e-2;
maxerr = abs(approxlevel - truelevel);
 


end
