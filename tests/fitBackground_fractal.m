function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Exponential background fit
%======================================================

t = linspace(0,5,100);
d = 3;

k = 0.5;
d = 3;
bckg = exp(-(k*t).^(d/3));
k = 1;
d = 2;
bckg2 = exp(-(k*t).^(d/3));
k = 1.5;
d = 4;
bckg3 = exp(-(k*t).^(d/3));

data2fit = bckg(1:end);
data2fit2 = bckg2(1:end);
data2fit3 = bckg3(1:end);

tfit = t(1:end);

[fit,results] = fitbackground(data2fit,t,tfit,'fractal');
[fit2,results2] = fitbackground(data2fit2,t,tfit,'fractal');
[fit3,results3] = fitbackground(data2fit3,t,tfit,'fractal');


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
  title(sprintf('k = %.3f (0.5) d = %.3f (3)',results.DecayRate,results.FractalDimension))
  subplot(132)
  plot(t,bckg2,t,fit2)
  title(sprintf('k = %.3f (1) d = %.3f (2)',results2.DecayRate,results2.FractalDimension))
  subplot(133)
  plot(t,bckg3,t,fit3)
  title(sprintf('k = %.3f (1.5) d = %.3f (4)',results3.DecayRate,results3.FractalDimension))
end

end