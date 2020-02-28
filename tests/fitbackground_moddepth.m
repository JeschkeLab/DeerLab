function [pass,maxerr] = test(opt)

% Check that fitbackground() can fit/accept modulation depth and use logfits

t = linspace(0,10,100);
lam0 = 0.25;
k = 0.05;
V = dipolarsignal(t,3,'Background',td_exp(t,k),'moddepth',lam0);

tstart = 2;

[~ ,lamfit] = fitbackground(V,t,@td_exp,tstart);
pass(1) = all(abs(lamfit - lam0) < 1e-2);

[~ ,~,kfit] = fitbackground(V,t,@td_exp,tstart,'ModDepth',lam0);
pass(2) = all(abs(kfit - k) < 1e-2);

[~ ,lamfit] = fitbackground(V,t,@td_exp,tstart,'Logfit',true);
pass(3) = all(abs(lamfit - lam0) < 1e-2);

[~ ,~,kfit] = fitbackground(V,t,@td_exp,tstart,'ModDepth',lam0,'Logfit',true);
pass(4) = all(abs(kfit - k) < 1e-2);

pass = all(pass);
maxerr = max(abs(lamfit - lam0));


end