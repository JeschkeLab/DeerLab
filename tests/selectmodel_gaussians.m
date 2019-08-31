function [err,data,maxerr] = test(opt,olddata)

%Test if selectmethod can identify that the optimal method is a two
%gaussian model as given as the input signal

Dimension = 200;
dt = 0.016;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.3 5 0.3 0.5];
Distribution = rd_twogaussian(r,InputParam);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

Models = {@rd_onegaussian,@rd_twogaussian,@rd_threegaussian};

[optimum,metric] = selectmodel(Models,DipEvoFcn,r,K,'aicc');

err = optimum~=2;
data = [];
maxerr = [];


if opt.Display
figure(8),clf
plot(metric)
end

end

