function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Exponential background fit
%======================================================
clear fitbackground

t = linspace(0,5,1000);
d = 3;
k = 0.5;
bckg = exp(-(k*t).^(d/3));
data2fit = bckg(1:end);

tstart = t(20);

tic
prefit = fitbackground(data2fit,t,@td_exp,tstart);
pre1 = toc;
tic
prefit2 = fitbackground(data2fit,t,@td_strexp,tstart);
pre2 = toc;
tic
postfit = fitbackground(data2fit,t,@td_exp,tstart);
post1 = toc;
tic
postfit2 = fitbackground(data2fit,t,@td_strexp,tstart);
post2 = toc;

err(1) = any(abs(prefit - postfit)>1e-15);
err(2) = any(abs(prefit2 - postfit2)>1e-15);
err(3) = post1 > pre1/4;
err(4) = post2 > pre2/4;

err = any(err);
maxerr = max(abs(prefit - postfit));
data = [];


end