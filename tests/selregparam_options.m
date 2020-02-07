function [err,data,maxerr] = test(opt,olddata)

N = 60;
dt = 0.008;
t = linspace(0,dt*N,N);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
S = K*P;

alphaopt1 = selregparam(S.',K,r,'tikhonov',{'aic','gml','gcv'},'NoiseLevel',1e-5);

alphaopt2 = selregparam(S,K,r,'tikhonov',{'aic','gml','gcv'},'RegOrder',2);
alphaopt3 = selregparam(S,K,r,'tikhonov',{'aic','gml','gcv'},'TolFun',1e-10);

alphaopt4 = selregparam(S,K,r,'huber',{'aic','gml','gcv'});
alphaopt5 = selregparam(S,K,r,'huber',{'aic','gml','gcv'},'HuberParameter',1.35);

alphaopt6 = selregparam(S,K,r,'huber',{'aic','gml','gcv'},'GlobalWeights',1);


err(1) = any(abs(alphaopt1 - alphaopt2) > 1e-2);
err(2) = any(abs(alphaopt1 - alphaopt3) > 1e-2);
err(3) = any(abs(alphaopt1 - alphaopt3) > 1e-2);
err(4) = any(abs(alphaopt4 - alphaopt5) > 1e-2);
err(5) = any(abs(alphaopt4 - alphaopt6) > 1e-2);

err = any(err);
data = [];
maxerr = 0;

end