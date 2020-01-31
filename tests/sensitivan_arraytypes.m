function [err,data,maxerr] = test(opt,olddata)

n = 20;
x = linspace(0,1,n);

param.a = linspace(1,4,2);
param.b = linspace(0.02,0.1,2);

stats  = sensitivan(@(param)myfcn(param,x),param);

err(1) = any(size(stats(1).mean)>1);
err(2) = any(size(stats(2).mean)~=[n 1]);
err(3) = any(size(stats(3).mean)~=[n n]);


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