function [err,data,maxerr] = test(opt,olddata)

%Test if selectmethod can identify that the optimal method is a two
%gaussian model as given as the input signal

Dimension = 500;
TimeStep = 0.016;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [3 0.3 5 0.3 0.5];
Distribution = InputParam(5)*gaussian(DistanceAxis,InputParam(1),InputParam(2))/(1/sqrt(2*pi)*1/InputParam(2)) ...
    + (1-InputParam(5))*gaussian(DistanceAxis,InputParam(3),InputParam(4))/(1/sqrt(2*pi)*1/InputParam(4));
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

Models = {@onegaussian,@twogaussian,@threegaussian,@onerice,@tworice};

tic
[preoptimum] = selectmodel(Models,DipEvoFcn,DistanceAxis,Kernel,'aicc');
pre  = toc;
tic
[postoptimum] = selectmodel(Models,DipEvoFcn,DistanceAxis,Kernel,'aicc');
post  = toc;

err(1) = any(postoptimum~=preoptimum);
err(2) = post > pre/4;
data = [];
maxerr = [];


if opt.Display
figure(8),clf
plot(metric)
end

end

