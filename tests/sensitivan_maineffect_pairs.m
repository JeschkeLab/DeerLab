function [pass,maxerr] = test(opt)

% Check that sensitivan computed the main effects correctly

x = linspace(0,100,100);

param.a = linspace(1,4,2);
param.b = linspace(0.02,0.1,2);

[~,fctrs]  = sensitivan(@(param)myfcn(param,x),param);
mainEffect = fctrs.main;

% Pass: main effects are correctly computed
pass = mainEffect.a < mainEffect.b;
pass = all(pass);
 
maxerr = NaN;

end


function y = myfcn(p,x)

a = p.a;
b = p.b;

rng(a)
y = whitegaussnoise(x,b);

end