function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Exponential background fit
%======================================================
clear fitbackground

t = linspace(0,5,1000);
d = 3;
k = 0.5;
bckg = exp(-(k*t).^(d/3));
data2fit = bckg(20:end);

tfit = t(20:end);

tic
prefit = fitbackground(data2fit,t,tfit,'exponential');
pre1 = toc;
tic
prefit2 = fitbackground(data2fit,t,tfit,'polyexp');
pre2 = toc;
tic
postfit = fitbackground(data2fit,t,tfit,'exponential');
post1 = toc;
tic
postfit2 = fitbackground(data2fit,t,tfit,'polyexp');
post2 = toc;

err(1) = any(abs(prefit - postfit)>1e-15);
err(2) = any(abs(prefit2 - postfit2)>1e-15);
err(3) = post1 > pre1/4;
err(4) = post2 > pre2/4;

err = any(err);
maxerr = max(abs(prefit - postfit));
data = [];


end