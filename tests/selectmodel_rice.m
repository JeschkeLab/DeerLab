function [err,data,maxerr] = test(opt,olddata)

%Test if selectmethod can identify that the optimal method is a two
%rice model as given as the input signal
Dimension = 300;
TimeStep = 0.016;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [3 0.2 5.5 0.3 0.5];
Distribution = rd_tworice(DistanceAxis,InputParam);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

Models = {@rd_onerice,@rd_tworice,@rd_threerice};

[optimum,metric] = selectmodel(Models,DipEvoFcn,DistanceAxis,Kernel,'aicc');

err = optimum~=2;
data = [];
maxerr = [];


if opt.Display
figure(8),clf
plot(metric)
end

end

