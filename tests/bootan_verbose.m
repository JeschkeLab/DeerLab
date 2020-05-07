function [pass,maxerr] = test(opt)

% Check that bootan's verbose mode works without crashes

sig = 0.05;

t = linspace(0,5,100);
r = linspace(2,6,100);
P = dd_gauss(r,[4 0.8]);
K = dipolarkernel(t,r);
V = K*P;

parfit = fitparamodel(V,@dd_gauss,r,K);

Vfit = K*dd_gauss(r,parfit) + whitegaussnoise(t,sig);
evalc('bootan(@bootfcn,V,Vfit,10,''verbose'',true)');

% Pass 1-2: function with verbose does not crash
pass = true;

pass = all(pass);

maxerr = NaN;

    function [parfit,Pfit] = bootfcn(V)
        parfit = fitparamodel(V,@dd_gauss,r,K);
        Pfit = dd_gauss(r,parfit);
    end


end