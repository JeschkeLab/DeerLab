function [pass,maxerr] = test(opt)

% Check whether regoperator returns correct matrices

n = 20;

L0 = regoperator(n,0);
L0ref = eye(n);
% Pass 1: the 0th order operator is equal to reference
pass(1) = isequal(L0,L0ref);

L1 = regoperator(n,1);
L1ref = diff(eye(n),1);
% Pass 2: the 1st order operator is equal to reference
pass(2) = isequal(L1,L1ref);

L2 = regoperator(n,2);
L2ref = diff(eye(n),2);
% Pass 3: the 2nd order operator is equal to reference
pass(3) = isequal(L2,L2ref);

L3 = regoperator(n,3);
L3ref = diff(eye(n),3);
% Pass 4: the 3rd order operator is equal to reference
pass(4) = isequal(L3,L3ref);

pass = all(pass);

maxerr = NaN;

end