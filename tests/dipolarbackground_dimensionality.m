function [pass,maxerr] = test(opt)

% Check indifference of dipolarbackground() towards input dimensionality

t = linspace(-1,5,150);
kappa = 0.3;
lam = 0.5;
Bmodel = @(t,lam) bg_exp(t,kappa,lam);
path(1,:) = [1-lam NaN];
path(2,:) = [lam 0];

B1 = dipolarbackground(t,path,Bmodel);
B2 = dipolarbackground(t.',path,Bmodel);


% Pass 1: all backgrounds are equal
pass(1) = isequal(B1,B2);
% Pass 2: all backgrounds are column vectors
pass(2) = iscolumn(B1) & iscolumn(B2);

pass = all(pass);

maxerr = NaN;
 

end