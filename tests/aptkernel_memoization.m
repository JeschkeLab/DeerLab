function [pass,maxerr] = test(opt)

% Check that aptkernel() memoization works

t = linspace(0,3,800);
% First run: slow
tic 
preK = aptkernel(t);
pre = toc;
% Second run: fast (cached)
tic 
postK = aptkernel(t);
post = toc;

% Pass 1: cached results should run at least 10x faster
err(1) = post < pre/10;
% Pass 2: cached results should be equal
err(2) = isequal(preK.Base,postK.Base);

pass = all(err);
 
maxerr = NaN;


end