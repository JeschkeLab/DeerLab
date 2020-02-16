function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check Tikhonov regularization
%=======================================
rng(1)

t = linspace(0,4,200);
r = time2dist(t);
P = rd_onegaussian(r,[4,0.5]);

K = dipolarkernel(t,r);
S = K*P;
noise = whitegaussnoise(t,0.01);
S = S + noise;

Pfit = fitregmodel(S,K,r,'tikhonov','aic');
err(1) = any(abs(Pfit - P)>8e-2);
maxerr = max(abs(Pfit - P));

err = any(err);
data = [];

if opt.Display
 	figure(8),clf
    subplot(121)
    hold on
    plot(t,S)
    plot(t,K*Pfit)
    subplot(122)
    hold on
    plot(r,P,'k') 
    plot(r,Pfit,'r')
end

end