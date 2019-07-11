function [err,data] = test(opt,olddata)

%=======================================
% Check TV regularization
%=======================================

Dimension = 200;
TimeStep = 0.008;
rmin = (4*TimeStep*52.04/0.85)^(1/3);
rmax = 6*(Dimension*TimeStep/2)^(1/3);
TimeAxis = linspace(0,TimeStep*Dimension,Dimension);
DistanceAxis = linspace(rmin,rmax,Dimension);
Distribution = gaussian(DistanceAxis,3,0.15);
Distribution = Distribution/sum(Distribution);

Kernel = getKernel(Dimension,TimeStep*1000);
Background = exp(-0.5*TimeAxis);

DipEvoFcn = Kernel*Distribution;

%Set optimal regularization parameter (found numerically lambda=0.13)
options = DAoptions('RegParam',0.005,'Solver','cvx','RegMatrixOrder',3);
TVResult1 = regularize(DipEvoFcn,Kernel,'tv',options);
err(1) = any(abs(TVResult1 - Distribution)>1e-3);
err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(DistanceAxis,Distribution,'k') 
    plot(DistanceAxis,TVResult1,'r')
end

end