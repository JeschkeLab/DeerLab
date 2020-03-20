function [pass,maxerr] = test(opt)

% Check indifference of selectmodel() towards input dimensionality

t = linspace(0,5,80);
S = dipolarsignal(t,3);
r = linspace(1,6,50);
K = dipolarkernel(t,r);

[opt1,f1,params1] = selectmodel({@dd_onegauss,@dd_fourgauss},S,r,K,'aicc');
[opt2,f2,params2] = selectmodel({@dd_onegauss,@dd_fourgauss},S.',r,K,'aicc');
[opt3,f3,params3] = selectmodel({@dd_onegauss,@dd_fourgauss},S,r.',K,'aicc');
[opt4,f4,params4] = selectmodel({@dd_onegauss,@dd_fourgauss},S.',r.',K,'aicc');

% Pass 1: all functionals are equal 
pass(1) = iscolumn(f1) & iscolumn(f2) & iscolumn(f3) & iscolumn(f4);
% Pass 2: all parameter vectors are rows 
pass(2) = ~iscolumn(params1) & ~iscolumn(params2) & ~iscolumn(params3) & ~iscolumn(params4);
% Pass 3: all optimas are equal 
pass(3) = isequal(opt1,opt2,opt3,opt4);
% Pass 4: all parameter vectors are equal 
pass(4) = isequal(params1,params2,params3,params4);

pass = all(pass);

maxerr = NaN;
 

end