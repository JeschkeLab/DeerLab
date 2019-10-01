function [err,data] = test(opt,olddata)

N = 20;
r = linspace(1,5,N);

L1 = regoperator(N,2);
L2 = regoperator(r,2);
  
err = any(abs(L1 - L2)>1e-10);
maxerr = max(abs(L1 - L2));
data = [];

end