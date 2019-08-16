function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check TV regularization
%=======================================
Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));
Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;


rESEEM = 1.5;
EseemFreq = 52.04/(rESEEM^3);
ESEEM = 0.5*(2 + exp(-7*TimeAxis).*cos(2*pi*EseemFreq*TimeAxis));
Signal = DipEvoFcn.*ESEEM';

Filtered = longpass(TimeAxis,Signal,1.5);

L = regoperator(Dimension,2);
RegParam = regparamrange(Kernel,L);
RegParam2 = selregparam(RegParam,Filtered,Kernel,L,'gml');
Result = fitregmodel(Filtered,Kernel,DistanceAxis,L,'tikhonov',RegParam2,'Solver','fnnls');

error = abs(Result - Distribution);
err(1) = any(error>3e-1);
maxerr = max(error);
err = any(err);
data = [];

if opt.Display
    figure(9),clf
    subplot(122),hold on
    plot(DistanceAxis,Distribution)
    plot(DistanceAxis,Result)
    Result = fitregmodel(Signal,Kernel,DistanceAxis,L,'tikhonov',RegParam2,'Solver','fnnls');
    plot(DistanceAxis,Result)
    subplot(121),hold on
    plot(TimeAxis,DipEvoFcn)
    plot(TimeAxis,Filtered)
    plot(TimeAxis,Signal)
end

end