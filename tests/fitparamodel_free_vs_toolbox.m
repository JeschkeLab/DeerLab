function [pass,maxerr] = test(opt)

% Check that fitparamodel() finds the same solution with a free solver and a toolbox solver

rng(1)

t = linspace(0,2,400);
r = time2dist(t);
par0  = [3 0.5];
P = dd_gauss(r,par0);
S = dipolarsignal(t,r,P,'noiselevel',0.05);
K = dipolarkernel(t,r);

param0 = [2 0.1];
[parfit1,Pfit1] = fitparamodel(S,@dd_gauss,r,K,param0,'Solver','fminsearchcon');
[parfit2,Pfit2] = fitparamodel(S,@dd_gauss,r,K,param0,'Solver','fmincon');

[parfit3,Pfit3] = fitparamodel(S,@dd_gauss,r,K,param0,'Solver','lsqnonlin');
[parfit4,Pfit4] = fitparamodel(S,@dd_gauss,r,K,param0,'Solver','fminsearchcon');

% Pass 1-2: fmincon (toolbox) and fminsearchcon (free) find the same solution
pass(1) = all(abs(Pfit1 - Pfit2) < 1e-3);
pass(2) = all(abs(parfit1 - parfit2) < 1e-3);

% Pass 3-4: lsqnonlin (toolbox) and fminsearchcon (free) find the same solution
pass(3) = all(abs(Pfit3 - Pfit4) < 1e-3);
pass(4) = all(abs(parfit3 - parfit4) < 1e-3);

pass = all(pass);

maxerr = max(abs(Pfit1 - Pfit2));
 
if opt.Display
   plot(r,P,'k',r,Pfit1,r,Pfit2,r,Pfit3,r,Pfit4)
   legend('truth','fit 1','fit 2','fit 3','fit 4')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end