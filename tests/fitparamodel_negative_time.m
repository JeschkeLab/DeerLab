function [err,data,maxerr] = test(opt,olddata)

t = linspace(-5,5,300);
r = time2dist(t);
P = rd_onegaussian(r,[4 0.4]);
K = dipolarkernel(t,r);
V = K*P;
V = V + whitenoise(300,0.02);
V = V/max(V);


Pfit = fitparamodel(V,K,r,@rd_onegaussian);

error = abs(Pfit - P);
err = any(error>7e-2);
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