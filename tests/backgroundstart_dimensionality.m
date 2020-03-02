function [pass,maxerr] = test(opt)

% Check indifference of backgroundstart() towards input dimensionality

t = linspace(0,5,30);
bckg = exp(-(0.5*t));
% Clear aptkernel since it is used internally
clear aptkernel
t1 = backgroundstart(bckg,t,@td_exp);
clear aptkernel
t2 = backgroundstart(bckg.',t,@td_exp);
clear aptkernel
t3 = backgroundstart(bckg,t.',@td_exp);
clear aptkernel
t4 = backgroundstart(bckg.',t.',@td_exp);

% Pass 1: all start times are equal
pass(1) = isequal(t1,t2,t3,t4);
% Pass 1: all values are scalars
pass(2) = isscalar(t1) & isscalar(t2) & isscalar(t3) & isscalar(t4);

pass = all(pass);

maxerr = NaN;
 

end