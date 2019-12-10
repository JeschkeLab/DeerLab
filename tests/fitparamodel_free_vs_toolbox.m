function [err,data,maxerr] = test(opt,oldata)

rng(1)

Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
paramin  = [3,0.5];
P = rd_onegaussian(r,paramin);

S = dipolarsignal(t,r,P,'noiselevel',0.05);

K = dipolarkernel(t,r);


param0 = [2 0.1];
[fitparam1,Pfit1] = fitparamodel(S,@rd_onegaussian,r,K,param0,'Solver','fminsearchcon');
[fitparam2,Pfit2] = fitparamodel(S,@rd_onegaussian,r,K,param0,'Solver','fmincon');

[fitparam3,Pfit3] = fitparamodel(S,@rd_onegaussian,r,K,param0,'Solver','lsqnonlin');
[fitparam4,Pfit4] = fitparamodel(S,@rd_onegaussian,r,K,param0,'Solver','nlsqbnd');

err(1) = any(abs(Pfit1 - Pfit2)>1e-5);
err(2) = any(abs(fitparam1 - fitparam2)>1e-3);

err(3) = any(abs(Pfit3 - Pfit4)>1e-5);
err(4) = any(abs(fitparam3 - fitparam4)>1e-3);

err = any(err);

maxerr = max(abs(Pfit1 - Pfit2));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,S,'k.')
   plot(t,K*Pfit1,'r')
   plot(t,K*Pfit2,'b')

end

end