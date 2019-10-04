function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.3]);

K = dipolarkernel(t,r);
S = dipolarsignal(t,r,P,'noiselevel',0.15);

alpha = 50000;

if opt.Display
    figure(8),clf
    axhandle = axes();
else
    axhandle = [];
end

Pfit = obir(S,K,r,'tikh',alpha,'NoiseLevelAim',0.05,'DivergenceStop',true,'Axishandle',axhandle);

err = any(abs(Pfit - P)>1e-1);
maxerr = max(abs(Pfit - P));
data = [];

if opt.Display
 	figure(8),clf
    hold on
    plot(r,P,'k') 
    plot(r,Pfit,'b')
    legend('truth','OBIR','Tikh')
end

end