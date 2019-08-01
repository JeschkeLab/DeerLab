

Dimension = 300;
TimeStep = 0.016;
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = time2dist(TimeAxis);
InputParam = [3 0.4 5.5 0.4 0.5];
Distribution = twogaussian(DistanceAxis,InputParam);
Distribution = Distribution/sum(Distribution);

Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

Models = {@onerice,@tworice,@threerice,@onegaussian,@twogaussian,@threegaussian};

[optimum,aicc,bic] = selectmodel(Models,DipEvoFcn,DistanceAxis,Kernel);

err = optimum~=2;
data = [];
maxerr = [];

for i=1:length(Models)
tags{i} = func2str(Models{i});
end

figure(8),clf
subplot(121)
plot(DistanceAxis,Distribution,'k','LineWidth',1.5)
grid on,box on
xlabel('Distance [nm]')
ylabel('P(r)')
subplot(122)
hold on
plot(aicc,'b-o','LineWidth',1.5)
plot(bic,'r-o','LineWidth',1.5)
legend('AICc','BIC')
grid on,box on,axis tight
ylabel('Metric')
set(gca,'xtick',1:length(Models),'xticklabel',tags)
xtickangle(45)
