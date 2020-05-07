function [pass,maxerr] = test(opt)

% Check indifference of selregparam() towards input dimensionality

t = linspace(-1,4,20);
r = linspace(1,6,50);
P = dd_gauss(r,[3 0.3]);
V = dipolarsignal(t,r,P,'noiselevel',0.02);
K = dipolarkernel(t,r);

[alpha1,F1,alphas1,res1,pen1] = selregparam(V,K,'tikh','aic');
[alpha2,F2,alphas2,res2,pen2] = selregparam(V.',K,'tikh','aic');


% Pass 1: all regularization parameter values are equal
pass(1) = isequal(alpha1,alpha2);
% Pass 2: all selection functionals are equal
pass(2) = isequal(F1,F2);
% Pass 3: all evaluated alpha values are equal
pass(3) = isequal(alphas1,alphas2);
% Pass 4: all evaluated alpha vectors are columns (if not empty)
if ~isempty(alphas1)
    pass(4) = iscolumn(alphas1) & iscolumn(alphas2);
else
    pass(4) = true;
end
% Pass 5: all evaluated functional vectors are columns
pass(5) = iscolumn(F1) & iscolumn(F2);
% Pass 6: all residual vectors are columns
pass(6) = iscolumn(res1) & iscolumn(res2);
% Pass 7: all penalty vectors are columns
pass(7) = iscolumn(pen1) & iscolumn(pen2);

pass = all(pass);

maxerr = NaN;
 

end