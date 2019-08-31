function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
Distribution = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;


%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.02;
RegMatrix = regoperator(Dimension,3);
RegFunctional = @(Dist)(1/2*norm(K*Dist - DipEvoFcn)^2 + RegParam^2*max(RegMatrix*Dist)^2);
Result = fitregmodel(DipEvoFcn,K,r,RegMatrix,RegFunctional,RegParam,'Solver','fmincon');

err = any(abs(Result - Distribution)>1e-1);
maxerr = max(abs(Result - Distribution));

data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,Distribution,'k') 
    plot(r,Result,'r')
end

end