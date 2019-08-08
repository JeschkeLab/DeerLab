function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Polynomial background fit
%======================================================

t = linspace(0,5,100);
d = 3;
k = 0.5;
bckg = polyval([1 1],t);
bckg2 = polyval([2 1 1],t);
bckg3 = polyval([3 0 -1 1],t);

data2fit = bckg(20:end);
data2fit2 = bckg2(20:end);
data2fit3 = bckg3(20:end);
tfit = t(20:end);

polyOrder = 1;
polyOrder2 = 2;
polyOrder3 = 3;

fit = fitbackground(data2fit,t,tfit,'polynomial',polyOrder);
fit2 = fitbackground(data2fit2,t,tfit,'polynomial',polyOrder2);
fit3 = fitbackground(data2fit3,t,tfit,'polynomial',polyOrder3);


err(1) = any(abs(fit - bckg)>1e-5);
err(2) = any(abs(fit2 - bckg2)>1e-5);
err(3) =  any(abs(fit3 - bckg3)>1e-5);
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
  legend('truth','fit')
end

end