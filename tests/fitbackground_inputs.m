function [pass,maxerr] = test(opt)


% Check that fitbackground() works with all input schemes

t = linspace(0,10,100);
lam0 = 0.25;
k = 0.05;
V = dipolarsignal(t,3,'Background',td_exp(t,k),'moddepth',lam0);
tstart = 4.2424;
tend = 10;

B1 = fitbackground(V,t,@td_exp,'ModDepth',0.25);
B2 = fitbackground(V,t,@td_exp,tstart,'ModDepth',0.25);
B3 = fitbackground(V,t,@td_exp,[tstart tend],'ModDepth',0.25);
B4 = fitbackground(V,t,@td_exp);

% Pass 1-3: all background are well fitted
pass(1) = all(abs(B1 - B2) < 1e-4);
pass(2) = all(abs(B3 - B2) < 1e-4);
pass(3) = all(abs(B4 - B2) < 1e-4);

pass = all(pass);

maxerr = max(abs(B1 - B2));
 

end