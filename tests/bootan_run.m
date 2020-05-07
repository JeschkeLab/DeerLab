function [pass,maxerr] = test(opt)

% Check that bootan() internal loop works

sig = 0.05;

t = linspace(0,5,200);
r = linspace(2,6,300);
P = dd_gauss(r,[4 0.8]);
K = dipolarkernel(t,r);
V = K*P;

parfit = fitparamodel(V,@dd_gauss,r,K);

Vfit = K*dd_gauss(r,parfit) + whitegaussnoise(t,sig);

bootan(@bootfcn,V,Vfit,10,'verbose',false);

% Pass 1: no crashes happened
pass = true;
 
maxerr = NaN;

    function parfit = bootfcn(V)
        parfit = fitparamodel(V,@dd_gauss,r,K);
    end


end