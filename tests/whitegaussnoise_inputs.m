function [pass,maxerr] = test(opt)

% Check the different input schemes for whitegaussnoise()

t = linspace(0,5,200);
rng(2)
noise1 = whitegaussnoise(200,0.05);
rng(2)
noise2 = whitegaussnoise(t,0.05);

% Pass: both input schemes lead to the same output 
pass = all(noise1 == noise2);

maxerr = NaN;
 


end