function [pass,maxerr] = test(opt)

Dimension = 100;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;

[OptParam,Functionals,RegParams] = selregparam(DipEvoFcn,K,r,'huber',{'aic','gcv'},'NonNegConstrained',false);

%Accept testif all values are the same (should be as there is no noise)
pass = all(any(OptParam - OptParam' > 1e-2));
maxerr = max(max(OptParam - OptParam'));
 

if opt.Display
   figure(8),clf
   hold on
   plot(RegParams,Functionals{1})
   plot(RegParams,Functionals{2})
end


end