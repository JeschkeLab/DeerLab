function [err,data,maxerr] = test(opt,olddata)

x = linspace(0,100,100);

param.a = linspace(1,4,2);
param.b = linspace(0.02,0.1,2);
param.c = linspace(0.04,0.3,2);

[~,~,~,mainEffect]  = sensitivan(@(param)myfcn(param,x),param);

err(1) = mainEffect{1}(1)>mainEffect{1}(2);
err(2) = mainEffect{2}(1)>mainEffect{2}(3);

err = any(err);
data = [];
maxerr = 0;

end


function [y,z] = myfcn(p,x)

a = p.a;
b = p.b;
c = p.c;

rng(a)
y = whitegaussnoise(x,b);
z = whitegaussnoise(x,c);

end