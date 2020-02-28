function [pass,maxerr] = test(opt)

t = linspace(0,5,200);

nu1 = time2freq(t);
nu2 = time2freq(t,200);

err(1) = any(abs(nu1-nu2)>1e-10);
pass = all(err);
maxerr = max(abs(nu1-nu2));
 

end