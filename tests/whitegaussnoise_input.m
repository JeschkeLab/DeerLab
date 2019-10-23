function [err,data,maxerr] = test(opt,olddata)

t = linspace(0,5,200);
rng(2)
N1 = whitegaussnoise(200,0.05);
rng(2)
N2 = whitegaussnoise(t,0.05);

err = any(N1~=N2);
maxerr = max(N1-N2);
data = [];


end