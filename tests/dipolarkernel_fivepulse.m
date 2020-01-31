
function [err,data,maxerr] = test(opt,olddata)


r = linspace(2,6,500);
P = rd_onegaussian(r,[3,0.8]);

tau1 = 4.24;
tau2 = 4.92;

N = 800;
t1 = linspace(0,10,N);
t2 = 0.3;

t = (tau1 + tau2) - (t1 + t2);
taus = [tau1 tau2];
ts = {t1 t2};

prob = 0.8;
V = td_dmpdeer(t,r,P,taus,ts,prob);
V = V/max(V);
t = fliplr(t);
V = flipud(V);

eta1 = 0;
eta2 = tau2 - t2;

lam1 = prob^2;
lam2 = prob*(1-prob);
lam0 = (1-prob)^2 + prob*(1-prob)*dipolarsignal(tau2-t2,r,P);
lambdas = [lam0 lam1 lam2];
etas = [eta1 eta2];

K = dipolarkernel(t,r,'multipathway',{lambdas etas});

err = any(abs(V - K*P)>1e-3);
maxerr = max(max(abs(V - K*P)));
data = [];

if opt.Display
   plot(t,V,t,K*P) 
end

end