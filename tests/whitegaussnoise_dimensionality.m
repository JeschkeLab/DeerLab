function [pass,maxerr] = test(opt)

% Check indifference of time2freq() towards input dimensionality

t = linspace(0,5,80);

rng(2)
noise1 = whitegaussnoise(t,0.05);
rng(2)
noise2 = whitegaussnoise(t.',0.05);

% Pass 1:  both noise vectors are column vectors
pass(1) = iscolumn(noise1) & iscolumn(noise2);
% Pass 2: both noise vectors are equal
pass(2) = isequal(noise1,noise2);

pass = all(pass);

maxerr = max(noise2 - noise1);
 

end