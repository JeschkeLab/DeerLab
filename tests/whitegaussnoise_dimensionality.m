function [pass,maxerr] = test(opt)

t = linspace(0,5,80);

rng(2)
noise1 = whitegaussnoise(t,0.05);
rng(2)
noise2 = whitegaussnoise(t.',0.05);

err(1) = ~iscolumn(noise1) | ~iscolumn(noise2);
err(2) = ~isequal(noise1,noise2);

pass = all(err);
maxerr = max(noise2 - noise1);
 

end