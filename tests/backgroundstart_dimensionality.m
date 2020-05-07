function [pass,maxerr] = test(opt)

% Check indifference of backgroundstart() towards input dimensionality

t = linspace(0,5,10);
bckg = exp(-(0.2*t));
% Clear aptkernel since it is used internally
t1 = backgroundstart(bckg,t,@bg_exp);
t2 = backgroundstart(bckg.',t,@bg_exp);
t3 = backgroundstart(bckg,t.',@bg_exp);
t4 = backgroundstart(bckg.',t.',@bg_exp);

% Pass 1: all start times are equal
pass(1) = isequal(t1,t2,t3,t4);
% Pass 1: all values are scalars
pass(2) = isscalar(t1) & isscalar(t2) & isscalar(t3) & isscalar(t4);

pass = all(pass);

maxerr = NaN;
 

end