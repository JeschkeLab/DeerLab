function [err,data,maxerr] = test(opt,olddata)

t = linspace(0,10,100);
lam0 = 0.25;
k = 0.05;
V = dipolarsignal(t,3,'Background',td_exp(t,k),'moddepth',lam0);

tstart = 2;

[~ ,lamfit] = fitbackground(V,t,@td_exp,tstart);
err(1) = any(abs(lamfit - lam0)>1e-2);

[~ ,~,kfit] = fitbackground(V,t,@td_exp,tstart,'ModDepth',lam0);
err(2) = any(abs(kfit - k)>1e-2);

[~ ,lamfit] = fitbackground(V,t,@td_exp,tstart,'Logfit',true);
err(3) = any(abs(lamfit - lam0)>1e-2);

[~ ,~,kfit] = fitbackground(V,t,@td_exp,tstart,'ModDepth',lam0,'Logfit',true);
err(4) = any(abs(kfit - k)>1e-2);

err = any(err);
maxerr = max(abs(lamfit - lam0));
data = [];

if opt.Display
  figure,clf
  hold on
  plot(t,td_exp(t,k),t,Bfit)
end

end