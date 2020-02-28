function [pass,maxerr] = test(opt)

% Check indifference of correctzerotime() towards input dimensionality

t = linspace(-1,4,100);
S = dipolarsignal(t,3);
t = t + abs(min(t));

t1 = correctzerotime(S,t);
t2 = correctzerotime(S.',t);
t3 = correctzerotime(S,t.');
t4 = correctzerotime(S.',t.');

% Pass 1: all corrected time axes are equal
pass(1) = isequal(t1,t2,t3,t4);
% Pass 2: all corrected time axes are column vectors
pass(2) = iscolumn(t1) & iscolumn(t2) & iscolumn(t3) & iscolumn(t4);

pass = all(pass);

maxerr = max(abs(t1 - t2));
 

end