function [pass,maxerr] = test(opt)

% Check the different input schemes of regoperator()

N = 20;
r = linspace(1,5,N);

L1 = regoperator(N,2);
L2 = regoperator(r,2);

% Pass: both schemes yield the same results
pass = isequal(L1,L2);

maxerr = NaN;
 

end