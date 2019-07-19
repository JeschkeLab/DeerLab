function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check TV regularization
%=======================================
Dimension = 200;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);
Kernel = getKernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;


rESEEM = 1.5;
EseemFreq = 52.04/(rESEEM^3);
ESEEM = 0.5*(2 + exp(-7*TimeAxis).*cos(2*pi*EseemFreq*TimeAxis));
Signal = DipEvoFcn.*ESEEM';

Filtered = longpassFilter(TimeAxis,Signal,2);

L = getRegMatrix(Dimension,2);
RegParam = getRegParamRange(Kernel,L);
RegParam2 = selectRegParam(RegParam,Filtered,Kernel,L,'gml');
Result = regularize(Filtered,Kernel,L,'tikhonov',RegParam2,'Solver','fnnls');

error = abs(Result - Distribution);
err(1) = any(error>1e-2);
maxerr = max(error);
err = any(err);
data = [];

if opt.Display
    figure(9),clf
    subplot(122),hold on
    plot(DistanceAxis,Distribution)
    plot(DistanceAxis,Result)
    subplot(121),hold on
    plot(TimeAxis,DipEvoFcn)
    plot(TimeAxis,Filtered)
    plot(TimeAxis,Signal)
end

end