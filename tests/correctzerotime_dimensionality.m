function [err,data,maxerr] = test(opt,olddata)


t = linspace(-1,4,100);
S = dipolarsignal(t,3);
t = t + abs(min(t));

t1 = correctzerotime(S,t);
t2 = correctzerotime(S.',t);
t3 = correctzerotime(S,t.');
t4 = correctzerotime(S.',t.');


err(1) = ~isequal(t1,t2,t3,t4);
err(2) = ~iscolumn(t1) | ~iscolumn(t2) | ~iscolumn(t3) |~iscolumn(t4);

err = any(err);
maxerr = max(abs(t1 - t2));
data = [];

end