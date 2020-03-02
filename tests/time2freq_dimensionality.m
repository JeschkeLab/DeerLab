function [pass,maxerr] = test(opt)

% Check indifference of time2freq() towards input dimensionality

t = linspace(0,5,80);

nu1 = time2freq(t);
nu2 = time2freq(t.');

% Pass 1:  both frequency axes are column vectors
pass(1) = iscolumn(nu1) & iscolumn(nu2);
% Pass 2: both frequency axes are equal
pass(2) = isequal(nu1,nu2);

pass = all(pass);

maxerr = NaN;
 

end