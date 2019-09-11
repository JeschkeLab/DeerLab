function [err,data,maxerr] = test(opt,olddata)

t = linspace(-3,5,200);
r = time2dist(t);
P = rd_onegaussian(r,[4 0.4]);
K = dipolarkernel(t,r);
V = K*P;
rng(2)
V = V + whitegaussnoise(200,0.01);
rng('default')
V = V/max(V);
alpha = 2;
Pfit = fitregmodel(V,K,r,'tikhonov',alpha);

error = abs(Pfit - P);
err = any(error>5e-1);
maxerr= max(error);
data = [];

if opt.Display
figure(8)
clf
subplot(121)
plot(r,P,r,Pfit)
subplot(122)
plot(t,V,t,K*Pfit)
end





end