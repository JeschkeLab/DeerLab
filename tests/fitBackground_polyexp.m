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

data2fit = bckg(20:end);
data2fit2 = bckg2(20:end);
data2fit3 = bckg3(20:end);

tfit = t(20:end);

fit = fitbackground(data2fit,t,tfit,'polyexp');
fit2 = fitbackground(data2fit2,t,tfit,'polyexp');
fit3 = fitbackground(data2fit3,t,tfit,'polyexp');

err(1) = any(abs(fit - bckg)>1e-5);
err(2) = any(abs(fit2 - bckg2)>1e-5);
err(3) =  any(abs(fit3 - bckg3)>1e-5);
err = any(err);
maxerr = max(fit - bckg);
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