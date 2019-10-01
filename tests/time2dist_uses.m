function [err,data,maxerr] = test(opt,oldata)

t = linspace(0,5,200);

r1 = time2dist(t);
r2 = time2dist(t,200);
[~,rmin,rmax] = time2dist(t);
r3 = linspace(rmin,rmax,200);

err(1) = any(abs(r1-r2)>1e-10);
err(2) = any(abs(r1-r3)>1e-10);
err = any(err);
maxerr = max(abs(r1-r2));
data = [];

end