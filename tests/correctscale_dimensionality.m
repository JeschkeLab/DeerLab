function [err,data,maxerr] = test(opt,olddata)


t = linspace(-1,4,100);
S = dipolarsignal(t,3,'scale',1e4);
t = t + abs(min(t));

V1 = correctscale(S,t);
V2 = correctscale(S.',t);
V3 = correctscale(S,t.');
V4 = correctscale(S.',t.');


err(1) = ~isequal(V1,V2,V3,V4);
err(2) = ~iscolumn(V1) | ~iscolumn(V2) | ~iscolumn(V3) |~iscolumn(V4);

err = any(err);
maxerr = max(abs(V1 - V2));
data = [];

end