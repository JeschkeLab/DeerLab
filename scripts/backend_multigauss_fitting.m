
Dimension = 300;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [4 0.2 4 1 3 0.4 0.4 0.4];
Distribution = threegaussian(DistanceAxis,InputParam);
Distribution = Distribution/(1/sqrt(2*pi)*1/InputParam(2));
Distribution = Distribution/sum(Distribution);
Noise = rand(Dimension,1);
Noise = Noise - mean(Noise);
Noise = 0.0*Noise/max(Noise);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution + Noise;

[FitDistribution,FitParam,optimum,metrics] = multigauss(DipEvoFcn,Kernel,DistanceAxis,8);

figure(8),clf
subplot(131)
plot(TimeAxis,DipEvoFcn,'k','LineWidth',1.5)
hold on
plot(TimeAxis,Kernel*FitDistribution,'r','LineWidth',1.5)
grid on,box on
legend('Truth','Fit')
xlabel('Time [\mus]')
ylabel('D(r)')
subplot(132)
plot(DistanceAxis,Distribution,'k','LineWidth',1.5)
hold on
plot(DistanceAxis,FitDistribution,'r','LineWidth',1.5)
grid on,box on
legend('Truth','Fit')
xlabel('Distance [nm]')
ylabel('P(r)')
subplot(133)
hold on
plot(metrics{1},'b-o','LineWidth',1.5)
plot(metrics{2},'r-o','LineWidth',1.5)
legend('AICc','BIC')
grid on,box on,axis tight
ylabel('Metric')
xlabel('N-Gaussian model')