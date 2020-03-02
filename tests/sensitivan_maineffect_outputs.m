function [pass,maxerr] = test(opt)

% Check that sensitivan computed the main effects correctly

x = linspace(0,100,100);

param.a = linspace(1,4,2);
param.b = linspace(0.02,0.1,2);
param.c = linspace(0.04,0.3,2);

[~,fctrs]  = sensitivan(@(param)myfcn(param,x),param);
mainEffect = fctrs.main;

% Pass: main effects are correctly computed
pass(1) = mainEffect{1}.a < mainEffect{1}.b;
pass(2) = mainEffect{2}.a < mainEffect{2}.c;

pass = all(pass);
 
maxerr = NaN;

end


function [y,z] = myfcn(p,x)

a = p.a;
b = p.b;
c = p.c;

rng(a)
y = whitegaussnoise(x,b);
z = whitegaussnoise(x,c);

end