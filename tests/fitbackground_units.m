function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Exponential background fit
%======================================================

t = linspace(0,5,100);
d = 3;

k = 0.5;
bckg = exp(-(k*t).^(d/3));
data2fit = bckg(1:end);

tstart = t(20);
%us
fit1 = fitbackground(data2fit,t,@td_exp,tstart);
%ns
t = t*1000;
fit2 = fitbackground(data2fit,t,@td_exp,tstart);

err = any(abs(fit1 - fit2)>1e-6);
err = any(err);
maxerr = max(abs(fit1 - fit2));
data = [];

if opt.Display
  figure,clf
  hold on
  plot(t,fit1)
  plot(t,fit2)
end

end