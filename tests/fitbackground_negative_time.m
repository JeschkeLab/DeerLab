function [err,data,maxerr] = test(opt,olddata)

t = linspace(-0.5,5,300);
r = time2dist(t);
P = rd_onegaussian(r,[4 0.4]);
K = dipolarkernel(t,r);
S = K*P;
rng(2)
S = S + whitegaussnoise(300,0.02);
rng('default')
S = S/max(S);
lambda = 0.2;
k = 0.3;
B = exp(-k*(abs(t)))';
F = (1-lambda) + lambda*S;
V = F.*B;
B = (1-lambda).*B;
tstart = backgroundstart(V,t,@td_exp);
[Bfit,lambdafit] = fitbackground(V,t,@td_exp,tstart);
Bfit = (1-lambdafit).*Bfit;

error = abs(B - Bfit);
err = any(error>1e-2);
maxerr= max(error);
data = [];

if opt.Display
figure(8)
clf
plot(t,V,t,B,t,Bfit)
end

end