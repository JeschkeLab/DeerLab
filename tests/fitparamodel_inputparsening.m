function [err,data,maxerr] = test(opt,oldata)

Dimension = 200;
t = linspace(0,5,Dimension);
r = time2dist(t);

P = rd_onegaussian(r,[3,0.5]);
B = td_strexp(t,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;

timemodel = @(t,param)K*rd_onegaussian(r,param);


InitialGuess = [2 0.1];
[~,Pfit1] = fitparamodel(S,@rd_onegaussian,r,K,InitialGuess,'Solver','lsqnonlin');
[~,Bfit1] = fitparamodel(B,@td_strexp,t,InitialGuess,'Solver','lsqnonlin');
[~,Bfit2] = fitparamodel(B,@td_strexp,t,'Solver','lsqnonlin');
fitparam = fitparamodel(S,timemodel,t,InitialGuess,'Solver','lsqnonlin');
[~,Pfit3] = fitparamodel(S,@rd_onegaussian,r,K,'Solver','lsqnonlin');

Pfit2 = rd_onegaussian(r,fitparam);

err(1) = any(abs(Pfit1 - P)>1e-5);
err(2) = any(abs(Bfit1 - B)>1e-5);
err(3) = any(abs(Bfit2 - B)>1e-5);
err(4) = any(abs(Pfit2 - P)>1e-5);
err(5) = any(abs(Pfit3 - P)>1e-5);

err = any(err);

maxerr = max(abs(Pfit1 - P));
data = [];

if opt.Display
   figure(1),clf,hold on
   plot(t,P,'b')
   plot(t,Pfit,'r')
end

end