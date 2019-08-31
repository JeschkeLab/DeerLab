function [err,data,maxerr] = test(opt,olddata)

%Test if selectmethod can identify that the optimal method is a two
%gaussian model as given as the input signal

Dimension = 500;
dt = 0.016;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.3 5 0.3 0.5];
Distribution = rd_twogaussian(r,InputParam);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

Models = {@rd_onegaussian,@rd_twogaussian,@rd_threegaussian,@rd_onerice,@rd_tworice};

tic
[preoptimum] = selectmodel(Models,DipEvoFcn,r,K,'aicc');
pre  = toc;
tic
[postoptimum] = selectmodel(Models,DipEvoFcn,r,K,'aicc');
post  = toc;

err(1) = any(postoptimum~=preoptimum);
err(2) = post > pre/3;
data = [];
maxerr = [];


end

