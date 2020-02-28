function [pass,maxerr] = test(opt)

t = linspace(0,5,200);
rng(2)
N1 = whitegaussnoise(200,0.05);
rng(2)
N2 = whitegaussnoise(t,0.05);

pass = all(N1~=N2);
maxerr = max(N1-N2);
 


end