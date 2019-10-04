function [err,data,maxerr] = test(opt,olddata)

dt = 0.5e-9;
nu1 = 0.5e9;
N = 300;
t = linspace(0,dt*N,N);
S = exp(-5e7*t).*(cos(2*pi*nu1*t));
sampl = 1/2*(1/dt);
wp = 0.01e9;
ws =  0.04e9;

[Sfilt1] = winlowpass(S,ws,wp,sampl);
[Sfilt2] = winlowpass(S,ws,wp,sampl,'MinimalAttenuation',50);


err = any(abs(Sfilt1 - Sfilt2)>1e-10);
maxerr = max(abs(Sfilt1 - Sfilt2));
data = [];

if opt.Display
figure(1),clf
subplot(121)
hold on
plot(t,S)
plot(t,Sfilt)

nu = time2freq(t,2*N);

subplot(122)
hold on
plot(nu,spec)
plot(nu,specfilt)
xline(wp,'k--');
end


end