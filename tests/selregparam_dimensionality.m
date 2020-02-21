function [err,data,maxerr] = test(opt,olddata)

t = linspace(-1,4,20);
r = linspace(1,6,50);
P = rd_onegaussian(r,[3 0.3]);
S = dipolarsignal(t,r,P,'noiselevel',0.02);
K = dipolarkernel(t,r);

[alpha1,F1,alphas1,res1,pen1] = selregparam(S,K,r,'tikh','aic');
[alpha2,F2,alphas2,res2,pen2] = selregparam(S.',K,r,'tikh','aic');
[alpha3,F3,alphas3,res3,pen3] = selregparam(S,K,r.','tikh','aic');
[alpha4,F4,alphas4,res4,pen4] = selregparam(S.',K,r.','tikh','aic');


err(1) = ~isequal(alpha1,alpha2,alpha3,alpha4);
err(2) = ~isequal(F1,F2,F3,F4);
err(3) = ~isequal(alphas1,alphas2,alphas3,alphas4);
err(4) = ~iscolumn(alphas1) | ~iscolumn(alphas2) | ~iscolumn(alphas3) | ~iscolumn(alphas4);
err(5) = ~iscolumn(F1) | ~iscolumn(F2) | ~iscolumn(F3) | ~iscolumn(F4);
err(6) = ~iscolumn(res1) | ~iscolumn(res2) | ~iscolumn(res3) | ~iscolumn(res4);
err(7) = ~iscolumn(pen1) | ~iscolumn(pen2) | ~iscolumn(pen3) | ~iscolumn(pen4);

err = any(err);

maxerr = max(abs(alpha1 - alpha2));
data = [];

end