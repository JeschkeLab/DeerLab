function [err,data,maxerr] = test(opt,olddata)


t = linspace(-1,4,100);

K1 = aptkernel(t);
clear aptkernel
K2 = aptkernel(t.');


err(1) = ~isequal(K1.Base,K2.Base);
err(2) = ~isequal(K1.FreqAxis,K2.FreqAxis);
err(3) = ~isequal(K1.t,K2.t);
err(4) = ~isequal(K1.Crosstalk,K2.Crosstalk);
err(5) = ~iscolumn(K1.t) | ~iscolumn(K2.t);
err(6) = ~iscolumn(K1.FreqAxis) | ~iscolumn(K2.FreqAxis);

err = any(err);
maxerr = max(abs(K1.FreqAxis - K2.FreqAxis));
data = [];

end