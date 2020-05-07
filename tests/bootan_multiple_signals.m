function [pass,maxerr] = test(opt)

% Check that bootan() internal loop works

sig = 0.05;

t1 = linspace(0,5,100);
t2 = linspace(0,2,150);
r = linspace(2,6,100);
P = dd_gauss(r,[4 0.8]);
K1 = dipolarkernel(t1,r);
V1 = K1*P + whitegaussnoise(t1,0.01);
K2 = dipolarkernel(t2,r);
V2 = K2*P + whitegaussnoise(t2,0.02);

parfit = fitparamodel({V1,V2},@dd_gauss,r,{K1,K2});
Vfit1 = K1*dd_gauss(r,parfit);
Vfit2 = K2*dd_gauss(r,parfit);

results = bootan(@bootfcn,{V1,V2},{Vfit1,Vfit2},10,'verbose',false);

% Pass 1-2: outputs has correct structure
pass(1) = numel(results)==2;
pass(2) = iscell(results);

pass = all(pass);

maxerr = NaN;

    function [parfit,Pfit] = bootfcn(Vs)
        parfit = fitparamodel(Vs,@dd_gauss,r,{K1,K2});
        Pfit = dd_gauss(r,parfit);
    end


end