function [err,data,maxerr] = test(opt,olddata)

x = linspace(0,50,50);

param.a = linspace(1,4,2);
param.b = linspace(0.02,0.1,2);

[mean]  = sensitivan(@(param)myfcn(param,x),param);

err(1) = any(size(mean{1})>1);
err(2) = any(size(mean{2})~=[1 50]);
err(3) = any(size(mean{3})~=[50 50]);


err = any(err);
data = [];
maxerr = 0;

end


function [x,y,z] = myfcn(p,in)

a = p.a;
b = p.b;

rng(a)

%scalar
x = a;
%vector
y = whitegaussnoise(in,b);
%matrix
z = y.*y.';


end