function [pass,maxerr] = test(opt)

% Check that fitbackground() can fit stretched exponentials


t = linspace(0,5,100).';
lam = 0.3;

k = 0.5;
d = 3;
bckg = (1-lam)*exp(-lam*k*(t).^(d/3));
k = 1;
d = 2;
bckg2 = (1-lam)*exp(-lam*k*(t).^(d/3));
k = 1.5;
d = 4;
bckg3 = (1-lam)*exp(-lam*k*(t).^(d/3));
data2fit = bckg(1:end);
data2fit2 = bckg2(1:end);
data2fit3 = bckg3(1:end);

[fit,lamfit] = fitbackground(data2fit,t,@bg_strexp,[min(t) max(t)]);
[fit2,lamfit2] = fitbackground(data2fit2,t,@bg_strexp,[min(t) max(t)]);
[fit3,lamfit3] = fitbackground(data2fit3,t,@bg_strexp,[min(t) max(t)]);

fit = (1-lamfit)*fit;
fit2 = (1-lamfit2)*fit2;
fit3 = (1-lamfit3)*fit3;

% Pass 1-3: all background are well fitted
pass(1) = all(abs(fit - bckg) < 1e-5);
pass(2) = all(abs(fit2 - bckg2) < 1e-5);
pass(3) = all(abs(fit3 - bckg3) < 1e-5);

pass = all(pass);
maxerr = max(abs(fit3 - bckg3));
 

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