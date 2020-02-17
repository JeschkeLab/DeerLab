function [err,data,maxerr] = test(opt,olddata)

t = linspace(0,5,80);

r1 = time2dist(t);
r2 = time2dist(t.');

err(1) = ~iscolumn(r1) | ~iscolumn(r2);
err(2) = ~isequal(r1,r2);

err = any(err);
maxerr = max(r1 - r2);
data = [];

end