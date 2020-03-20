function [pass,maxerr] = test(opt)

% Test global fit using Huber regularization

rng(2)
r = linspace(1,6,300);
P = dd_twogauss(r,[2,0.3,4,0.3,0.5]);

t1 = linspace(0,0.8,100);
K1 = dipolarkernel(t1,r);
S1 = K1*P + whitegaussnoise(t1,0.03);

t2 = linspace(0,1.6,200);
K2 = dipolarkernel(t2,r);
S2 = K2*P + whitegaussnoise(t2,0.05);

t3 = linspace(0,2.4,300);
K3 = dipolarkernel(t3,r);
S3 = K3*P + whitegaussnoise(t3,0.1);

regparam = 4;
Ss = {S1,S2,S3};
Ks = {K1,K2,K3};

Pglobal = fitregmodel(Ss,Ks,r,'huber',regparam,'Solver','fnnls');
Plocal1 = fitregmodel(S1,K1,r,'huber',regparam,'Solver','fnnls');
Plocal2 = fitregmodel(S2,K2,r,'huber',regparam,'Solver','fnnls');
Plocal3 = fitregmodel(S3,K3,r,'huber',regparam,'Solver','fnnls');

error_global = norm(P - Pglobal);
error_local1 = norm(P - Plocal1);
error_local2 = norm(P - Plocal2);
error_local3 = norm(P - Plocal3);

% Pass: the global fit is better than the local fits
pass = all(error_global < [error_local2 error_local3]);

maxerr = max(abs(Pglobal - P));

if opt.Display
    subplot(121)
    hold on
    plot(t1,S1,'k', t2,S2+1,'k', t3,S3+2,'k')
    plot(t1,K1*Pglobal,'r' ,t2,K2*Pglobal + 1,'r' ,t3,K3*Pglobal + 2,'r')
    xlabel('t [\mus]')
    ylabel('V(t)')
    legend('data','fit')
    grid on, axis tight, box on
    subplot(122)
    plot(r,P,'k', r,Pglobal,'c', r,Plocal1,'g--', r,Plocal2,'b--', r,Plocal3,'r--')
    xlabel('r [nm]')
    ylabel('P(r) [nm^{-1}]')
    legend('truth','global','local 1','local 2','local 3')
    grid on, axis tight, box on
end

end