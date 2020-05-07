function [pass,maxerr] = test(opt)

% Clear persistent variable
clear dipolarkernel

t = linspace(0,6,1000);
r = time2dist(t);

% First run: slow
tic
preK = dipolarkernel(t,r);
precached = toc;
% Second run: fast (cached)
tic
postK = dipolarkernel(t,r);
postcached = toc;

% Pass 1: cached results should run at least 10x faster
pass(1) = postcached <= precached/10;
% Pass 2: cached results should be equal
pass(2) = isequal(postK,preK);

pass = all(pass);

delta = abs(postK - preK);
maxerr = max(delta(:));

end