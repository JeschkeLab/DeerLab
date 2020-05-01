function [pass,maxerr] = test(opt)

% Check that fitbackground() can fit a sum of stretched exponentials


t = linspace(0,5,100).';
d = 1;
lam = 0.3;

k = 0.5;
k2 = 0.2;
bckg = 0.4*exp(-lam*k*(t).^(d)) + 0.6*exp(-lam*k2*(t).^(d));
bckg = (1-lam)*bckg;

data2fit = bckg(1:end);

tstart = t(20);

[fit,lamfit] = fitbackground(data2fit,t,@bg_sumstrexp,tstart);

fit = (1-lamfit)*fit;

% Pass 1-3: all background are well fitted
pass(1) = all(abs(fit - bckg) < 2e-3);

pass = all(pass);
maxerr = max(abs(fit - bckg));
 

%Plot results
if opt.Display
    plot(t,bckg,t,fit)
    legend('truth','fit')
    xlabel('t [\mus]')
    ylabel('B(t)')
    grid on, axis tight, box on
end

end