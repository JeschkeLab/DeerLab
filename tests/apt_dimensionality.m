function [err,data,maxerr] = test(opt,olddata)


t = linspace(-1,4,100);
S = dipolarsignal(t,3);
K = aptkernel(t);

Pfit1 = apt(S,K);
Pfit2 = apt(S.',K);


err(1) = ~isequal(Pfit1,Pfit2);
err(2) = ~iscolumn(Pfit1) | ~iscolumn(Pfit2);

err = any(err);
maxerr = max(abs(Pfit1 - Pfit2));
data = [];

end