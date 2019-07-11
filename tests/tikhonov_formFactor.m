function [err,data] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
TimeStep = 0.008;
rmin = (4*TimeStep*52.04/0.85)^(1/3);
rmax = 6*(Dimension*TimeStep/2)^(1/3);
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = linspace(rmin,rmax,Dimension);
Distribution = gaussian(DistanceAxis,3,0.5);
Distribution = Distribution/sum(Distribution);

Kernel = getKernel(Dimension,TimeStep*1000);
Background = exp(-0.15*TimeAxis)';

DipEvoFcn = Kernel*Distribution;
Cluster = (DipEvoFcn+5).*Background;
Background = Background./Cluster(1);

options = DAoptions('Background',Background);
Signal = DAsignal('TimeAxis',TimeAxis,'ExpData',Cluster);
Signal = Signal.prepare(options);

%Set optimal regularization parameter (found numerically lambda=0.13)
options = DAoptions('RegParam',0.13);

Kernel = getKernel(Dimension,TimeStep*1000,[],[],Background);


options = DAoptions('RegParam',0.13,'Solver','fnnls');
TikhResult1 = regularize(Signal.ClusterFcn,Kernel,'tikhonov',options);
options.Solver = 'bppnnls';
TikhResult2 = regularize(Signal.ClusterFcn,Kernel,'tikhonov',options);
options.Solver = 'lsqnonneg';
options.nonNegLSQsolTol = 1e-15;
TikhResult3 = regularize(Signal.ClusterFcn,Kernel,'tikhonov',options);

err(1) = any(abs(TikhResult1 - Distribution)>2e-6);
err(2) = any(abs(TikhResult2 - Distribution)>2e-6);
err(3) = any(abs(TikhResult3 - Distribution)>2e-4);
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,TikhResult,'r')
end

end