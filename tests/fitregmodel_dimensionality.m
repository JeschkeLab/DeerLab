function [pass,maxerr] = test(opt)



t = linspace(0,5,50);
r = linspace(1,6,50);
P = rd_onegaussian(r,[4 0.3]);

K = dipolarkernel(t,r);
S = K*P;

Pfit1 = fitregmodel(S,K,r);
Pfit2 = fitregmodel(S.',K,r);
Pfit3 = fitregmodel(S,K,r.');
Pfit4 = fitregmodel(S.',K,r.');

err(1) = ~isequal(Pfit1,Pfit2,Pfit3,Pfit4);
err(2) = ~iscolumn(Pfit1) | ~iscolumn(Pfit2) | ~iscolumn(Pfit3) | ~iscolumn(Pfit4);
pass = all(err);

maxerr = max(abs(Pfit1 - Pfit2));
 


end