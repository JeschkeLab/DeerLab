function [err,data,maxerr] = test(opt,olddata)

t = linspace(0,10,100);
lam0 = 0.25;
k = 0.05;
V = dipolarsignal(t,3,'Background',td_exp(t,k),'moddepth',lam0);

tstart = 4.2424;
tend = 10.0000;


B1 = fitbackground(V,t,@td_exp,[tstart tend],'Solver','lsqnonlin');
B2 = fitbackground(V,t,@td_exp,[tstart tend],'Solver','fminsearchcon');
B3 = fitbackground(V,t,@td_exp,[tstart tend],'Solver','nlsqbnd');

err(1) = any(abs(B1 - B2)>1e-6);
err(2) = any(abs(B3 - B2)>1e-6);

err = any(err);
maxerr = max(abs(B1 - B3));
data = [];

if opt.Display
  figure,clf
  hold on
  plot(t,B1,t,B2,t,B3)
end

end