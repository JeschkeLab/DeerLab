function [pass,maxerr] = test(opt)

% Check that fitbackground() can fit a product of stretched exponentials

t = linspace(0,5,100).';
d = 3;

k = 0.5;
k2 = 0.2;
bckg = exp(-(k*t).^(d/3)).*exp(-(k2*t).^(d/3));
k = 1;
k2 = 0.6;
bckg2 = exp(-(k*t).^(d/3)).*exp(-(k2*t).^(d/3));
k = 1.5;
k2 = 2;
bckg3 = exp(-(k*t).^(d/3)).*exp(-(k2*t).^(d/3));

data2fit = bckg(1:end);
data2fit2 = bckg2(1:end);
data2fit3 = bckg3(1:end);

tstart = t(20);

fit = fitbackground(data2fit,t,@bg_prodstrexp,tstart);
fit2 = fitbackground(data2fit2,t,@bg_prodstrexp,tstart);
fit3 = fitbackground(data2fit3,t,@bg_prodstrexp,tstart);

% Pass 1-3: all background are well fitted
pass(1) = all(abs(fit - bckg) < 1e-3);
pass(2) = all(abs(fit2 - bckg2) < 1e-3);
pass(3) =  all(abs(fit3 - bckg3) < 1e-3);

pass = all(pass);

maxerr = max(fit - bckg);
 

%Plot results
if opt.Display
    subplot(131)
    plot(t,bckg,t,fit)
    legend('truth','fit')
    xlabel('t [\mus]')
    ylabel('B(t)')
    grid on, axis tight, box on
    subplot(132)
    plot(t,bckg2,t,fit2)
    legend('truth','fit')
    xlabel('t [\mus]')
    ylabel('B(t)')
    grid on, axis tight, box on
    subplot(133)
    plot(t,bckg3,t,fit3)
    legend('truth','fit')
    xlabel('t [\mus]')
    ylabel('B(t)')
    grid on, axis tight, box on
end


end