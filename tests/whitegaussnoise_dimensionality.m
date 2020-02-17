function [err,data,maxerr] = test(opt,olddata)

t = linspace(0,5,80);

rng(2)
noise1 = whitegaussnoise(t,0.05);
rng(2)
noise2 = whitegaussnoise(t.',0.05);

err(1) = ~iscolumn(noise1) | ~iscolumn(noise2);
err(2) = ~isequal(noise1,noise2);

err = any(err);
maxerr = max(noise2 - noise1);
data = [];

end