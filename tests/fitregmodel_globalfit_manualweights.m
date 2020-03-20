function [pass,maxerr] = test(opt)


% Test global fit using Tikhonov regularization and user-given global weigths

rng(2)
r = linspace(1,6,300);
P = dd_twogaussian(r,[2,0.3,4,0.3,0.5]);

t1 = linspace(0,0.8,100);
K1 = dipolarkernel(t1,r);
S1 = K1*P + whitegaussnoise(t1,0.03);

t2 = linspace(0,1.6,200);
K2 = dipolarkernel(t2,r);
S2 = K2*P + whitegaussnoise(t2,0.05);

t3 = linspace(0,2.4,300);
K3 = dipolarkernel(t3,r);
S3 = K3*P + whitegaussnoise(t3,0.1);

regparam = 1;
Ss = {S1,S2,S3};
Ks = {K1,K2,K3};


Pglobal = fitregmodel(Ss,Ks,r,'tikh',regparam,'Solver','fnnls','GlobalWeights',[1 2 1]);

% Pass: the distribution is somewhat fitted, but it runs
pass = all(abs(P - Pglobal) < 5.5e-1);

maxerr = max(abs(P - Pglobal));

if opt.Display
    plot(r,P,'k', r,Pglobal,'r')
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    legend('truth','global','local 1','local 2','local 3')
    grid on, axis tight, box on
end

end