function [pass,maxerr] = test(opt)

t = linspace(0,5,80);
S = dipolarsignal(t,3);

X1 = winlowpass(S,3e6,2e6,8e6);
X2 = winlowpass(S.',3e6,2e6,8e6);

err(1) = ~iscolumn(X2) | ~iscolumn(X1);
err(2) = ~isequal(X2,X1);

pass = all(err);
maxerr = max(X2 - X1);
 

end