function [pass,maxerr] = test(opt)

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

approxlevel1 = noiselevel(S);
approxlevel2 = noiselevel(S,1/5*N);

err = abs(approxlevel1 - approxlevel2)>1e-10;
maxerr = abs(approxlevel1 - approxlevel2);
 


end
