function [err,data,maxerr] = test(opt,olddata)

%Test if selectmethod can identify that the optimal method is a two
%rice model as given as the input signal
Dimension = 300;
dt = 0.016;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
InputParam = [3 0.2 5.5 0.3 0.5];
Distribution = rd_tworice(r,InputParam);
Distribution = Distribution/sum(Distribution)/mean(diff(r));

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

Models = {@rd_onerice,@rd_tworice,@rd_threerice};

[optimum,metric] = selectmodel(Models,DipEvoFcn,r,K,'aicc');

err = optimum~=2;
data = [];
maxerr = [];


if opt.Display
figure(8),clf
plot(metric)
end

end

