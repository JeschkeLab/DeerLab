function [pass,maxerr] = test(opt)

% Check that obir() works with Tikhonov regularization using the fmincon solver

rng(1)
t = linspace(0,3,200);
r = linspace(0,5,100);
P = dd_gauss2(r,[2,0.3,3.5,0.3,0.5]);
K = dipolarkernel(t,r);
noiselvl = 0.1;
S = K*P + whitegaussnoise(t,noiselvl);

if opt.Display
    figure(8),clf
    axhandle = axes();
else
    axhandle = [];
end

alpha = 3;
Pobir = obir(S,K,r,'tikhonov',alpha,'NoiseLevelAim',noiselvl,'Solver','fmincon','axishandle',axhandle);

Preg = fitregmodel(S,K,r,'tikhonov',alpha);

pass = mean(abs(Pobir - P)) < mean(abs(Preg - P));

maxerr = max(abs(Pobir - P));
 
if opt.Display
    plot(r,P,'k',r,Preg,r,Pobir)
    legend('truth','regularization','OBIR')
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    grid on, axis tight, box on
end

end