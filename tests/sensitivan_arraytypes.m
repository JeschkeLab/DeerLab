function [pass,maxerr] = test(opt)

% Check that sensitivan works with different array types

n = 20;
x = linspace(0,1,n);

param.a = linspace(1,4,2);
param.b = linspace(0.02,0.1,2);

stats  = sensitivan(@(param)myfcn(param,x),param);

% Pass 1: first variable is a scalar
pass(1) = any(size(stats(1).mean) == 1);
% Pass 1: second variable is a vector
pass(2) = any(size(stats(2).mean) == [n 1]);
% Pass 3: third variable is a matrix
pass(3) = any(size(stats(3).mean) == [n n]);

pass = all(pass);
 
maxerr = NaN;

end


function [x,y,z] = myfcn(p,in)

a = p.a;
b = p.b;

rng(a)

%scalar
x = a;
%vector
y = whitegaussnoise(in,b);
%matrix
z = y.*y.';


end