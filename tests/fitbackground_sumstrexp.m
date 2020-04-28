function [pass,maxerr] = test(opt)

% Check that fitbackground() can fit a sum of stretched exponentials


t = linspace(0,5,100).';
d = 3;
lam = 0.3;

k = 0.5;
k2 = 0.2;
bckg = 0.4*exp(-k*(t).^(d/3)) + 0.6*exp(-k2*(t).^(d/3));
bckg = (1-lam)*bckg;
k = 1;
k2 = 0.6;
bckg2 = 0.4*exp(-k*(t).^(d/3)) + 0.6*exp(-k2*(t).^(d/3));
bckg2 = (1-lam)*bckg2;
k = 1.5;
k2 = 2;
bckg3 = 0.4*exp(-k*(t).^(d/3)) + 0.6*exp(-k2*(t).^(d/3));
bckg3 = (1-lam)*bckg3;

data2fit = bckg(1:end);
data2fit2 = bckg2(1:end);
data2fit3 = bckg3(1:end);

tstart = t(20);

[fit,lamfit] = fitbackground(data2fit,t,@bg_sumstrexp,tstart);
[fit2,lamfit2] = fitbackground(data2fit2,t,@bg_sumstrexp,tstart);
[fit3,lamfit3] = fitbackground(data2fit3,t,@bg_sumstrexp,tstart);

fit = (1-lamfit)*fit;
fit2 = (1-lamfit2)*fit2;
fit3 = (1-lamfit3)*fit3;

% Pass 1-3: all background are well fitted
pass(1) = all(abs(fit - bckg) < 2e-3);
pass(2) = all(abs(fit2 - bckg2) < 3e-3);
pass(3) = all(abs(fit3 - bckg3) < 2.5e-3);

pass = all(pass);
maxerr = max(abs(fit - bckg));
 

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