function [err,data,maxerr] = test(opt,olddata)

t = linspace(-5,1,400);
r = time2dist(t);
P = rd_onegaussian(r,[4,0.2]);
S = dipolarkernel(t,r)*P;
zt = abs(min(t));

[ct,czt] = correctzerotime(S,t+zt);


err(1) = any(abs(ct - t)>1e-10);
err(2) = abs(czt' - zt)>1e-10;

err = any(err);
data = [];
maxerr = max(abs(ct - t));

if opt.Display
   figure(8),clf
   hold on
   plot(ct,S,'.')
end

end