function [pass,maxerr] = test(opt)

% Check that fitbackground() can fit polynomials

clear fitbackground
t = linspace(0,3,100).';
bckg = polyval([-1 1],t);
bckg2 = polyval([-1 -1 1],t);
bckg3 = polyval([-1 -1 -1 1],t);

data2fit = bckg(1:end);
data2fit2 = bckg2(1:end);
data2fit3 = bckg3(1:end);
tstart = t(1);

[fit,lambda1] = fitbackground(data2fit,t,@bg_poly1,tstart);
[fit2,lambda2] = fitbackground(data2fit2,t,@bg_poly2,tstart);
[fit3,lambda3] = fitbackground(data2fit3,t,@bg_poly3,tstart);

fit = fit*(1-lambda1);
fit2 = fit2*(1-lambda2);
fit3 = fit3*(1-lambda3);

% Pass 1-3: all background are well fitted
pass(1) = all(abs(fit - bckg) < 1e-8);
pass(2) = all(abs(fit2 - bckg2) < 1e-8);
pass(3) =  all(abs(fit3 - bckg3) < 1e-8);

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