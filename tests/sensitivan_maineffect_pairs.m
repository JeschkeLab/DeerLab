function [err,data,maxerr] = test(opt,olddata)

x = linspace(0,100,100);

param.a = linspace(1,4,2);
param.b = linspace(0.02,0.1,2);

[~,~,~,mainEffect]  = sensitivan(@(param)myfcn(param,x),param);

err = mainEffect(1)>mainEffect(2);
err = any(err);
data = [];
maxerr = 0;

end


function y = myfcn(p,x)

a = p.a;
b = p.b;

rng(a)
y = whitegaussnoise(x,b);

end