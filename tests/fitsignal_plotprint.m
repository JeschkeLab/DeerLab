function [pass,maxerr] = test(opt)

% Check that fitsignal() plot&print works

rng(1)

t = linspace(0,5,100);
r = linspace(2,6,150);
P = dd_gauss(r,[4.5 0.6]);

info = ex_4pdeer(t);
parIn = [info.parameters.default];
pathinfo = ex_4pdeer(t,parIn);

kappa = 0.4;
Bmodel = @(t,lam) bg_exp(t,kappa,lam);
K = dipolarkernel(t,r,pathinfo,Bmodel);
V = K*P + whitegaussnoise(t,0.01);

if opt.Display
    set(0,'DefaultFigureVisible','on')
else
    set(0,'DefaultFigureVisible','off')
end

try
    evalc('fitsignal(V,t,r,''P'',@bg_exp,@ex_4pdeer)');
    pass = true;
catch
    pass = false;
end
set(0,'DefaultFigureVisible','on')

maxerr = NaN;


end