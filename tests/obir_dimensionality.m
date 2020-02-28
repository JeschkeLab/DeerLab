function [pass,maxerr] = test(opt)

rng(2)
t = linspace(-1,4,20);
r = linspace(1,6,30);
P = rd_onegaussian(r,[3 0.3]);
S = dipolarsignal(t,r,P,'noiselevel',0.02);
K = dipolarkernel(t,r);

Pfit1 = obir(S,K,r,'tikh','aic','MaxIter',5);
Pfit2 = obir(S.',K,r,'tikh','aic','MaxIter',5);
Pfit3 = obir(S,K,r.','tikh','aic','MaxIter',5);
Pfit4 = obir(S.',K,r.','tikh','aic','MaxIter',5);


err(1) = ~isequal(Pfit1,Pfit2,Pfit3,Pfit4);
err(2) = ~iscolumn(Pfit1) | ~iscolumn(Pfit2) | ~iscolumn(Pfit3) | ~iscolumn(Pfit4);

pass = all(err);

maxerr = max(abs(Pfit1 - Pfit2));
 

end