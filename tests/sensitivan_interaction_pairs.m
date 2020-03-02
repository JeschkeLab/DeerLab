function [pass,maxerr] = test(opt)

% Check that sensitivan computes interactions between factors correctly

x = linspace(0,100,100);

param.a = linspace(1,4,2);
param.b = linspace(0.02,0.1,2);
param.c = linspace(0.04,0.3,2);

[~,factors]  = sensitivan(@(param)myfcn(param,x),param);
Interaction = factors.inter;

% Pass 1-2: interaction values are correctly computed
pass(1) = Interaction{1}(1,3) < Interaction{1}(1,2);
pass(2) = Interaction{2}(1,2) < Interaction{2}(1,3);

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