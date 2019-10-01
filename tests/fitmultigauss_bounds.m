function [err,data,maxerr] = test(opt,oldata)


Dimension = 200;
t = linspace(0,3,200);
r = time2dist(t);
InputParam = [4 0.1 4.5 0.4 0.35];
P = rd_twogaussian(r,InputParam);

K = dipolarkernel(t,r);
S = K*P;

[Pfit1,param,nGaussOpt,metrics,Peval] = fitmultigauss(S,K,r,5,'aic','Lower',[2 0.05],'Upper',[6 0.5]);
[Pfit2] = fitmultigauss(S,K,r,5,'aic','Lower',[2 0.05]);
[Pfit3] = fitmultigauss(S,K,r,5,'aic','Upper',[6 0.5]);



err(1) = any(abs(Pfit1 - P)>7e-5);
err(2) = any(abs(Pfit1 - P)>7e-5);
err(3) = any(abs(Pfit1 - P)>7e-5);
err = any(err);

maxerr = max(abs(Pfit1 - P));
data = [];

if opt.Display
   figure(1),clf
   subplot(121)
   hold on
   plot(t,S,'b')
   plot(t,K*FitP,'r')
   subplot(122)
   hold on
   plot(r,P,'k',r,Pfit1,r,Pfit2,r,Pfit3)
end

end