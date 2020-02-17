function [err,data,maxerr] = test(opt,oldata)



t = linspace(0,5,100);
r = linspace(1,6,100);
P = rd_onegaussian(r,[4 0.3]);

K = dipolarkernel(t,r);
S = K*P;

[Pfit1,parfit1] = fitmultigauss(S,K,r,3,'aic');
[Pfit2,parfit2] = fitmultigauss(S.',K,r.',3,'aic');
[Pfit3,parfit3] = fitmultigauss(S.',K,r,3,'aic');
[Pfit4,parfit4] = fitmultigauss(S,K,r.',3,'aic');

err(1) = ~isequal(Pfit1,Pfit2,Pfit3,Pfit4);
err(2) = ~iscolumn(Pfit1) | ~iscolumn(Pfit2) | ~iscolumn(Pfit3) | ~iscolumn(Pfit4);
err(3) = iscolumn(parfit1) | iscolumn(parfit2) | iscolumn(parfit3) | iscolumn(parfit4);

err = any(err);

maxerr = max(abs(Pfit1 - Pfit2));
data = [];


end