clc
%Load experimental data
[ExpTimeAxis,ExpSignal] = eprload('CT_DEER_mix_28_36.DSC');

%Correct phase and zero-time
Signal = correctPhase(ExpSignal);
[Signal,TimeAxis] = correctZeroTime(Signal,ExpTimeAxis);
%Convert to us and truncate
TimeAxis = TimeAxis/1000;
[~,ZeroTimePos] = min(abs(TimeAxis));
Signal = Signal(ZeroTimePos:end);
TimeAxis = TimeAxis(ZeroTimePos:end);

%Fit background
FitStartPos = 20;
Data2fit = Signal(FitStartPos:end);
FitTimeAxis = TimeAxis(FitStartPos:end);
Background = fitBackground(Data2fit,TimeAxis,FitTimeAxis,'exponential');

%Ghost distance suppression
Signal = supressGhostDistances(Signal,2);

%Perform partial background correction
Signal = Signal./sqrt(Background);

%Prepare regularization
DistanceAxis = time2dist(TimeAxis);
Kernel = getKernel(TimeAxis,DistanceAxis,Background);
RegMatrix = getRegMatrix(length(Kernel),2);
RegParamRange = getRegParamRange(Kernel,RegMatrix);
RegParam = selectRegParam(RegParamRange,Signal,Kernel,RegMatrix,{'gml','lr'});

%Run regularization
Distribution1 = regularize(Signal,Kernel,RegMatrix,'tikhonov',RegParam(1),'Solver','fnnls');
Distribution2 = regularize(Signal,Kernel,RegMatrix,'tikhonov',RegParam(2),'Solver','fnnls');

figure(1),clf

subplot(121),hold on
plot(TimeAxis,Signal,'k')
plot(TimeAxis,sqrt(Background))
plot(TimeAxis,Kernel*Distribution1,'r')
plot(TimeAxis,Kernel*Distribution2,'b')
xlabel('Time [\mus]')
ylabel('Intensity [a.u.]')
legend('Exp','Bckg','GML','L-Curve')

subplot(122),hold on
plot(DistanceAxis,Distribution1,'r')
plot(DistanceAxis,Distribution2,'b')
xlabel('Distance [nm]')
ylabel('P(r)')
legend('GML','L-Curve')
