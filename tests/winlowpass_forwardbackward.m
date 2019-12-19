function [err,data,maxerr] = test(opt,olddata)

dt = 0.5e-9;
nu1 = 0.5e9;
N = 300;
t = linspace(0,dt*N,N);
signal = exp(-5e7*t).*(cos(2*pi*nu1*t));
% signal = signal/signal(1);
sampl = 1/2*(1/dt);
% ws  = 4e8;
wp = 0.01e9;
ws =  0.04e9;


filteredS = winlowpass(signal,ws,wp,sampl,'ForwardBackward',true);
filteredS = winlowpass(signal,ws,wp,sampl,'ForwardBackward',false);

spec = abs(fftshift(fft(signal,2*N)));
filteredspec = abs(fftshift(fft(filteredS,2*N)));

attenuation = mag2db(max(filteredspec)/max(spec));

err = attenuation>-30;
maxerr = attenuation;
data = [];

if opt.Display
figure(1),clf
subplot(121)
hold on
plot(t,signal)
plot(t,filteredS)

nu = time2freq(t,2*N);

subplot(122)
hold on
plot(nu,spec)
plot(nu,filteredspec)
xline(wp,'k--');
end


end