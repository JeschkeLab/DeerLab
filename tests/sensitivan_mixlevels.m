function [pass,maxerr] = test(opt)

% Check that sensitivan can deal with mixed-level factorial analysis

x = linspace(0,100,100);

param.a = linspace(1,4,2);
param.b = linspace(0.02,0.1,4);
param.c = linspace(0.04,0.3,3);

[~,fctrs]  = sensitivan(@(param)myfcn(param,x),param);
mainEffect = fctrs.main;

% Pass 1-2: main effects are correctly computed
pass(1) = mainEffect{1}.a < mainEffect{1}.b(1);
pass(2) = mainEffect{2}.a < mainEffect{2}.c(1);
% Pass 3: the output has the right dimension
pass(3) = numel(mainEffect{1}.c) == 2;

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