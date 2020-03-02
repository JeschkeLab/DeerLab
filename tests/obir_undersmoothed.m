function [pass,maxerr] = test(opt)

% Check that obir() enforces oversmoothing at the start

rng(1)
t = linspace(0,3,200);
r = linspace(0,5,100);
P = rd_twogaussian(r,[2,0.3,3.5,0.3,0.5]);
K = dipolarkernel(t,r);
noiselvl = 0.05;
S = K*P + whitegaussnoise(t,noiselvl);

alpha = 0.005;

if opt.Display
    figure(8),clf
    axhandle = axes();
else
    axhandle = [];
end

Pfit = obir(S,K,r,'tikh',alpha,'NoiseLevelAim',noiselvl,'Axishandle',axhandle);

% Pass: the distribution is well fitted
pass = all(abs(Pfit - P) < 7e-1);

maxerr = max(abs(Pfit - P));
 
 
if opt.Display
    plot(r,P,'k',r,Pfit)
    legend('truth','regularization','OBIR')
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    grid on, axis tight, box on
end

end