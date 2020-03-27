function [pass,maxerr] = test(opt)

% Check the different input schemes of fitparamodel()

t = linspace(0,5,200);
r = linspace(2,6,100);

P = dd_onegauss(r,[3,0.5]);
B = bg_strexp(t,[0.5,3]);
K = dipolarkernel(t,r);
S = K*P;

timemodel = @(t,param)K*dd_onegauss(r,param);
warning('off')

InitialGuess = [2 0.1];
[~,Pfit1] = fitparamodel(S,@dd_onegauss,r,K,InitialGuess,'Solver','lsqnonlin');
[~,Bfit1] = fitparamodel(B,@bg_strexp,t,InitialGuess,'Solver','lsqnonlin');
[~,Bfit2] = fitparamodel(B,@bg_strexp,t,'Solver','lsqnonlin');
fitparam = fitparamodel(S,timemodel,t,InitialGuess,'Solver','lsqnonlin');
[~,Pfit3] = fitparamodel(S,@dd_onegauss,r,K,'Solver','lsqnonlin');

warning('on')

Pfit2 = dd_onegauss(r,fitparam);

% Pass 1-5: all fits are consistent with the inputs
pass(1) = all(abs(Pfit1 - P) < 1e-5);
pass(2) = all(abs(Bfit1 - B) < 1e-5);
pass(3) = all(abs(Bfit2 - B) < 1e-5);
pass(4) = all(abs(Pfit2 - P) < 1e-5);
pass(5) = all(abs(Pfit3 - P) < 1e-5);

pass = all(pass);

maxerr = max(abs(Pfit1 - P));
 
if opt.Display
   plot(r,P,'k',r,Pfit1,'r')
   legend('truth','fit')
   xlabel('r [nm]')
   ylabel('P(r) [nm^{-1}]')
   grid on, axis tight, box on
end


end