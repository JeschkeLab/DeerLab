function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Polynomial background fit
%======================================================
clear fitbackground
t = linspace(0,3,100);
bckg = polyval([-1 1],t);
bckg2 = polyval([-1 -1 1],t);
bckg3 = polyval([-1 -1 -1 1],t);

data2fit = bckg(1:end);
data2fit2 = bckg2(1:end);
data2fit3 = bckg3(1:end);
tstart = t(1);

[fit,lambda1] = fitbackground(data2fit,t,@td_poly1,tstart);
[fit2,lambda2] = fitbackground(data2fit2,t,@td_poly2,tstart);
[fit3,lambda3] = fitbackground(data2fit3,t,@td_poly3,tstart);

fit = fit*(1-lambda1);
fit2 = fit2*(1-lambda2);
fit3 = fit3*(1-lambda3);

err(1) = any(abs(fit - bckg)>1e-8);
err(2) = any(abs(fit2 - bckg2)>1e-8);
err(3) =  any(abs(fit3 - bckg3)>1e-8);
err = any(err);
maxerr = max(abs(fit - bckg));
data = [];

if opt.Display
  figure(8),clf
  subplot(131)
  plot(t,bckg,t,fit)
  subplot(132)
  plot(t,bckg2,t,fit2)
  subplot(133)
  plot(t,bckg3,t,fit3)
  legend('truth','fit')
end

end