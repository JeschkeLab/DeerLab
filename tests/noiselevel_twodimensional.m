function [pass,maxerr] = test(opt)

% Test noiselevel() with 2D-datasets

rng(12);

t = linspace(0,5,300);
r = linspace(2,6,200);
P = dd_gauss(r,[4 0.4]);
K = dipolarkernel(t,r);

sigma0 = 0.1;
N = 500;
V = zeros(numel(t),N);
for i = 1:N
    V(:,i) = K*P + whitegaussnoise(t,sigma0);
end

sigma = noiselevel(V);

% Pass: noiselevels are equal
pass = abs(sigma - sigma0)<1e-3;

maxerr = max(abs(sigma - sigma0));
 

end