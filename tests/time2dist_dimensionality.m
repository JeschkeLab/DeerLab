function [pass,maxerr] = test(opt)

% Check indifference of time2dist() towards input dimensionality

t = linspace(0,5,80);

r1 = time2dist(t);
r2 = time2dist(t.');

% Pass 1:  both distance axes are column vectors
pass(1) = iscolumn(r1) & iscolumn(r2);
% Pass 2: both distance axes are equal
pass(2) = isequal(r1,r2);

pass = all(pass);

maxerr = NaN;
 

end