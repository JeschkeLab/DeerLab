function [pass,maxerr] = test(opt)

clear regularize

%=======================================
% Check TV regularization
%=======================================
Dimension = 80;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

%Set optimal regularization parameter (found numerically lambda=0.005)
RegParam = 0.001;
Result = fitregmodel(DipEvoFcn,K,r,'huber',RegParam,'Solver','fnnls','HuberParam',1.35);

error = abs(Result - P);
err(1) = any(error>1e-2);
maxerr = max(error);
pass = all(err);
 

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Result,'r')
end

end