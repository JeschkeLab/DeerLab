function [err,data] = test(opt,olddata)

%======================================================
% Zero-time correction function
%======================================================

t = linspace(-2,5,500);
r = time2dist(t);
P = rd_onegaussian(r,[4,0.2]);
S = dipolarkernel(t,r)*P;
zt = abs(min(t));

[ct,czt] = correctzerotime(S,t+zt);


err(1) = any(abs(ct - t)>1e-10);
err(2) = abs(czt' - zt)>1e-10;

err = any(err);
data = [];

if opt.Display
   figure(8),clf
   hold on
   plot(ct,S,'.')
end

end