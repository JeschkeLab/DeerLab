function [pass,maxerr] = test(opt)

% Check indifference of aptkernel() towards input dimensionality

t = linspace(-1,4,100);
K1 = aptkernel(t);
% Clear to remove cached results
clear aptkernel
K2 = aptkernel(t.');

% Pass 1: Kernel bases are equal
pass(1) = isequal(K1.Base,K2.Base);
% Pass 2: Frequency axes are equal
pass(2) = isequal(K1.FreqAxis,K2.FreqAxis);
% Pass 3: Time axes are equal
pass(3) = isequal(K1.t,K2.t);
% Pass 4: Crosstalk matrices are equal
pass(4) = isequal(K1.Crosstalk,K2.Crosstalk);
% Pass 5: Time axes are column vectors
pass(5) = iscolumn(K1.t) & iscolumn(K2.t);
% Pass 6: Frequency axes are column vectors
pass(6) = iscolumn(K1.FreqAxis) & iscolumn(K2.FreqAxis);

pass = all(pass);

maxerr = NaN;

end