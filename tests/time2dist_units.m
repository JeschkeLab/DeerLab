function [err,data,maxerr] = test(opt,oldata)

t = linspace(0,5,200);

%us
r1 = time2dist(t);
%ns
t = 1000*t;
r2 = time2dist(t);

err = any(abs(r1-r2)>1e-10);
maxerr = max(abs(r1-r2));
data = [];

end