function [pass,maxerr] = test(opt)

% Check whether regoperator returns correct matrices

r = linspace(2,5,100);
dr = mean(diff(r));

L0 = regoperator(r,0);
L0ref = eye(numel(r))/dr^0;
% Pass 1: the 0th order operator is equal to reference
pass(1) = isequal(L0,L0ref);

L1 = regoperator(r,1);
L1ref = diff(eye(numel(r)),1)/dr^1;
% Pass 2: the 1st order operator is equal to reference
pass(2) = isequal(L1,L1ref);

L2 = regoperator(r,2);
L2ref = diff(eye(numel(r)),2)/dr^2;
% Pass 3: the 2nd order operator is equal to reference
pass(3) = isequal(L2,L2ref);

L3 = regoperator(r,3);
L3ref = diff(eye(numel(r)),3)/dr^3;
% Pass 4: the 3rd order operator is equal to reference
pass(4) = isequal(L3,L3ref);

pass = all(pass);

maxerr = NaN;

end