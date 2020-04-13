function [pass,maxerr] = test(opt)

% Check that bootan() internal loop works

sig = 0.05;

t = linspace(0,5,100);
r = linspace(2,6,100);
P = dd_gauss(r,[4 0.8]);
K = dipolarkernel(t,r);
V = K*P;

parfit = fitparamodel(V,@dd_gauss,r,K);

Vfit = K*dd_gauss(r,parfit) + whitegaussnoise(t,sig);

results = bootan(@bootfcn,V,Vfit,'samples',80,'verbose',false);

% Pass 1-2: outputs has correct structure
pass(1) = numel(results)==2;
pass(2) = iscell(results);

pass = all(pass);

maxerr = NaN;

    function [parfit,Pfit] = bootfcn(V)
        parfit = fitparamodel(V,@dd_gauss,r,K);
        Pfit = dd_gauss(r,parfit);
    end


end