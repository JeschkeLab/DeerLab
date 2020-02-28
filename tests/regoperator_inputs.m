function [pass,maxerr] = test(opt)

N = 20;
r = linspace(1,5,N);

L1 = regoperator(N,2);
L2 = regoperator(r,2);
  
pass = all(abs(L1 - L2)>1e-10);
maxerr = max(max(abs(L1 - L2)));
 

end