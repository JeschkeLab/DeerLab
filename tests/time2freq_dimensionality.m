function [pass,maxerr] = test(opt)

t = linspace(0,5,80);

nu1 = time2freq(t);
nu2 = time2freq(t.');

err(1) = ~iscolumn(nu1) | ~iscolumn(nu2);
err(2) = ~isequal(nu1,nu2);

pass = all(err);
maxerr = max(nu1 - nu2);
 

end