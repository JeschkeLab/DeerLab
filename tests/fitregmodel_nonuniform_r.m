function [pass,maxerr] = test(opt)

% Check that fitregmodel works with non-uniform distributions

rng(1)

t = linspace(0,4,400);
r = sqrt(linspace(1,7^2,200));

P = dd_gauss(r,[4.5 0.5]);

K = dipolarkernel(t,r);
sigma = 0.01;
V = K*P + whitegaussnoise(t,sigma);

[Pfit] = fitregmodel(V,K,r,'tikh','aic');
Vfit = K*Pfit;

% Pass 1: non-uniform distribution is well fitted
pass(1) = all(abs(P - Pfit) < 1e-1);
% Pass 1: normalization works allright
pass(2) = round(Vfit(1),4) == 1;

pass = all(pass);
maxerr = max(P - Pfit);

if opt.Display
   plot(r,P,'ok',r,Pfit,'.-r')
   axis tight,grid on, box on
   xlabel('r [nm]'),ylabel('P(r) [nm^{-1}]')
end

end