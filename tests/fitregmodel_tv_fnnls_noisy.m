function [pass,maxerr] = test(opt)

%=======================================
% Check Tikhonov regularization
%=======================================

rng(2)
Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;
Noise = whitegaussnoise(Dimension,0.02);
%Set optimal regularization parameter (found numerically lambda=0.13)
RegParam = 0.1;
Resultfnnls = fitregmodel(DipEvoFcn+Noise,K,r,'tv',RegParam,'Solver','fnnls','RegOrder',3);

pass = all(abs(Resultfnnls - P)>2e-1);

maxerr = max(abs(Resultfnnls - P));

 

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k')
    plot(r,Resultfnnls,'b') 
    axis tight
end

end