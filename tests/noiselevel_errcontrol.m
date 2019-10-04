function [err,data,maxerr] = test(opt,olddata)

N = 200;
dt = 0.008;
t = linspace(0,dt*N,N);
r = time2dist(t);
InputParam = [3 0.5];
P = rd_onegaussian(r,InputParam);
P = P/sum(P);

K = dipolarkernel(t,r);
S = K*P;

rng(2)
Noise = rand(N,1);
Noise = Noise - mean(Noise);
Noise = 0.02*Noise/Noise(1);

S = S + Noise;

try
    noiselevel(S,1000);
    err = true;
catch
   err = false; 
end

maxerr = 0;
data = [];


end
