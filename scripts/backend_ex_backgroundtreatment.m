clc
%Load experimental data
[ExpTimeAxis,ExpSignal] = eprload('CT_DEER_mix_28_36.DSC');
Noise = rand(400,1);
Noise = Noise - mean(Noise);
Noise = 0.01*Noise/max(Noise)';
ExpSignal = ExpSignal/max(ExpSignal);
ExpSignal = ExpSignal + Noise;

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
Background = fitBackground(Data2fit,TimeAxis,FitTimeAxis);

%Ghost distance suppression
Signal = supressGhostDistances(Signal);

%Perform partial background correction
Signal2Reg = Signal./sqrt(Background);

%Prepare regularization
DistanceAxis = time2dist(TimeAxis);
Kernel = getKernel(TimeAxis,DistanceAxis,Background);
RegMatrix = getRegMatrix(length(Kernel),2);
RegParamRange = getRegParamRange(Kernel,RegMatrix);
RegParam = selectRegParam(RegParamRange,Signal2Reg,Kernel,RegMatrix,{'gml'});

%Run regularization
Distribution1 = regularize(Signal2Reg,Kernel,RegMatrix,'tikhonov',RegParam(1),'Solver','fnnls');

Signal2Reg2 = Signal;

%Prepare regularization
Kernel = getKernel(TimeAxis,DistanceAxis,Background,'KernelBType','full');
RegParamRange = getRegParamRange(Kernel,RegMatrix);
RegParam = selectRegParam(RegParamRange,Signal2Reg2,Kernel,RegMatrix,{'gml'});

Distribution2 = regularize(Signal2Reg2,Kernel,RegMatrix,'tikhonov',RegParam(1),'Solver','fnnls');

ModDepth = 1/Background(1) - 1;
Signal2Reg3 = Signal./Background;
Signal2Reg3 = (Signal2Reg3 - (1-ModDepth))/ModDepth - 1;
% Signal2Reg3 = Signal2Reg3/Signal2Reg3(1);
%Prepare regularization
Kernel = getKernel(TimeAxis,DistanceAxis);
RegParamRange = getRegParamRange(Kernel,RegMatrix);
RegParam = selectRegParam(RegParamRange,Signal2Reg3,Kernel,RegMatrix,{'gml'});

Distribution3 = regularize(Signal2Reg3,Kernel,RegMatrix,'tikhonov',RegParam(1),'Solver','fnnls');


figure(1),clf

subplot(121),hold on
plot(TimeAxis,Signal2Reg3,'k')
plot(TimeAxis,Kernel*Distribution1)
plot(TimeAxis,Kernel*Distribution2)
plot(TimeAxis,Kernel*Distribution3)

xlabel('Time [\mus]')
ylabel('Intensity [a.u.]')
legend('Sqrt(B)','B','division')

subplot(122),hold on
plot(DistanceAxis,Distribution1)
plot(DistanceAxis,Distribution2)
plot(DistanceAxis,Distribution3)

xlabel('Distance [nm]')
ylabel('P(r)')
legend('Sqrt(B)','B','division')
