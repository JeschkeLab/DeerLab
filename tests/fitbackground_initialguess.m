function [pass,maxerr] = test(opt)

% Check that fitbackground() works with passed initial guess values


t = linspace(0,5,100).';
d = 3;

k = 0.5;
bckg = exp(-(k*t).^(d/3));
k = 1;
bckg2 = exp(-(k*t).^(d/3));
k = 1.5;
bckg3 = exp(-(k*t).^(d/3));

data2fit = bckg(1:end);
data2fit2 = bckg2(1:end);
data2fit3 = bckg3(1:end);
tstart = t(20);

param1 = [0.45,2.8];
param2 = [0.8,3];
param3 = [1.1,3];

fit = fitbackground(data2fit,t,@bg_exp,tstart,'initialguess',param1);
fit2 = fitbackground(data2fit2,t,@bg_exp,tstart,'initialguess',param2);
fit3 = fitbackground(data2fit3,t,@bg_exp,tstart,'initialguess',param3);

% Pass 1-3: all background are well fitted
err(1) = any(abs(fit - bckg) < 1e-5);
err(2) = any(abs(fit2 - bckg2) < 1e-5);
err(3) = any(abs(fit3 - bckg3) < 1e-5);

pass = all(err);

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