function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Exponential background fit
%======================================================

t = linspace(0,5,100);
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

fit = fitbackground(data2fit,t,@td_exp,tstart);
fit2 = fitbackground(data2fit2,t,@td_exp,tstart);
fit3 = fitbackground(data2fit3,t,@td_exp,tstart);

err(1) = any(abs(fit - bckg)>1e-5);
err(2) = any(abs(fit2 - bckg2)>1e-5);
err(3) = any(abs(fit3 - bckg3)>1e-5);
err = any(err);
maxerr = max(abs(fit - bckg));
data = [];

if opt.Display
  figure,clf
  subplot(131)
  plot(t,bckg,t,fit)
  subplot(132)
  plot(t,bckg2,t,fit2)
  subplot(133)
  plot(t,bckg3,t,fit3)
end

end