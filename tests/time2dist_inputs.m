function [pass,maxerr] = test(opt)

% Check the different input schemes for time2dist()

t = linspace(0,5,200).';

r1 = time2dist(t);
r2 = time2dist(t,200);

% Pass: both input schemes lead to the same output 
pass = all(abs(r1-r2) < 1e-10);

maxerr = NaN;
 

end