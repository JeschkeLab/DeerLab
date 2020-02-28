function [pass,maxerr] = test(opt)

N = 200;
t = linspace(0,5,N);
r = time2dist(t);
K = dipolarkernel(t,r);
rng(2)
L = rand(1,N);
try
    regparamrange(K,L);
    err = false;
catch
    err = true;
end

maxerr = 0;
 

end