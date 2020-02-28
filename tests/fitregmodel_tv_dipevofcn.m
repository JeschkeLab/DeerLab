function [pass,maxerr] = test(opt)

%=======================================
% Check TV regularization
%=======================================
Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

%Set optimal regularization parameter (found numerically lambda=0.005)
RegParam = 0.0005;
TVResult1 = fitregmodel(DipEvoFcn,K,r,'tv',RegParam,'Solver','fmincon','RegOrder',3);

error = abs(TVResult1 - P);
err(1) = any(error>2e-1);
maxerr = max(error);
pass = all(err);
 

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,TVResult1,'r')
end

end