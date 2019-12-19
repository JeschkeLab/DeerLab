function [err,data,maxerr] = test(opt,olddata)

%Test if selectmethod can identify that the optimal method is a two
%gaussian model as given as the input signal

Dimension = 100;
dt = 0.016;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.3 5 0.3 0.5];
P = rd_twogaussian(r,InputParam);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

Models = {@rd_onegaussian,@rd_twogaussian,@rd_threegaussian};

[optimum,metric] = selectmodel(Models,DipEvoFcn,r,K,{'aic','aicc','bic'});
[optimum,metric] = selectmodel(Models,DipEvoFcn,r,K,'all');

err(1) = optimum(1)~=optimum(3);
err(2) = optimum(1)~=optimum(2);
err = any(err);
data = [];
maxerr = NaN;


if opt.Display
figure(8),clf
plot(metric)
end

end

