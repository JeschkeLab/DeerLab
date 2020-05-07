function [pass,maxerr] = test(opt)


% Check that fitbackground() works with all input schemes

t = linspace(0,5,100);
lam0 = 0.25;
k = 0.05;
B = bg_exp(t,k);
r = linspace(1,6,300);
P = dd_gauss(r,[3,0.3]);
V = dipolarsignal(t,r,P,lam0);
tstart = 2;
tend = 5;

B1 = fitbackground(V,t,@bg_exp,'ModDepth',lam0);
B2 = fitbackground(V,t,@bg_exp,tstart,'ModDepth',lam0);
B3 = fitbackground(V,t,@bg_exp,[tstart tend],'ModDepth',lam0);
B4 = fitbackground(V,t,@bg_exp);

% Pass 1-3: all background are well fitted
pass(1) = all(abs(B1 - B2) < 1e-4);
pass(2) = all(abs(B3 - B2) < 1e-4);
pass(3) = all(abs(B4 - B2) < 1e-4);

pass = all(pass);

maxerr = max(abs(B1 - B2));
 

end