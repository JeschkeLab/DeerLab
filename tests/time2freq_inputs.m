function [pass,maxerr] = test(opt)

% Check the different input schemes for time2freq()

t = linspace(0,5,200);

nu1 = time2freq(t);
nu2 = time2freq(t,200);

% Pass: both input schemes lead to the same output 
pass = all(abs(nu1-nu2) < 1e-10);

maxerr = NaN;
 

end