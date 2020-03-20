function [pass,maxerr] = test(opt)

% Check that fitparamodel() finds the same solution with a free solver and a toolbox solver

t = linspace(0,2,400);
r = time2dist(t);
paramin  = [3 0.5];
P = dd_onegaussian(r,paramin);
rng(1)
S = dipolarsignal(t,r,P,'noiselevel',0.05);
K = dipolarkernel(t,r);

param0 = [2 0.1];
[parafit1,Pfit1] = fitparamodel(S,@dd_onegaussian,r,K,param0,'Solver','fminsearchcon');
[parafit2,Pfit2] = fitparamodel(S,@dd_onegaussian,r,K,param0,'Solver','fmincon');

[fitparam3,Pfit3] = fitparamodel(S,@dd_onegaussian,r,K,param0,'Solver','lsqnonlin');
[fitparam4,Pfit4] = fitparamodel(S,@dd_onegaussian,r,K,param0,'Solver','fminsearchcon');

% Pass 1-2: fmincon (toolbox) and fminsearchcon (free) find the same solution
err(1) = all(abs(Pfit1 - Pfit2) < 1e-5);
err(2) = all(abs(parafit1 - parafit2) < 1e-3);
% Pass 3-4: lsqnonlin (toolbox) and fminsearchcon (free) find the same solution
err(3) = all(abs(Pfit3 - Pfit4) < 1e-5);
err(4) = all(abs(fitparam3 - fitparam4) < 1e-3);

pass = all(err);

maxerr = max(abs(Pfit1 - Pfit2));
 
if opt.Display
   plot(r,P,'k',r,Pfit1,r,Pfit2,r,Pfit3,r,Pfit4)
   legend('truth','fit 1','fit 2','fit 3','fit 4')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end

end