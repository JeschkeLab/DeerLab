function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;


%Set optimal regularization parameter (found numerically lambda=0.13)
tic
RegParam = 0.2;
RegMatrix = regoperator(Dimension,3);
RegFunctional = @(P)(1/2*norm(K*P - DipEvoFcn)^2 + RegParam^2*max(RegMatrix*P)^2);
Result = fitregmodel(DipEvoFcn,K,r,RegFunctional,RegParam,'Solver','fmincon');
toc

err = any(abs(Result - P)>5e-1);
maxerr = max(abs(Result - P));

data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Result,'r')
end

end