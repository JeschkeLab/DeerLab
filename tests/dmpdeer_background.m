function [err,data,maxerr] = test(opt,olddata)

r = linspace(2,6,500);
P = rd_onegaussian(r,[3,0.8]);

tau1 = 4.24;
tau2 = 4.92;

N = 800;
t1 = linspace(0,8,N);
t2 = 0.3;

t = (tau1 + tau2) - (t1 + t2);

taus = [tau1 tau2];
ts = {t1 t2};

prob = 0.8;

V1 = td_dmpdeer(t,r,P,taus,ts,prob);

prob = [0.8 0.8];
labels = [1 2];

Bparam = [0 3];

V2 = td_dmpdeer(t,r,P,taus,ts,prob,Bparam,labels);

err = any(abs(V1 - V2)>1e-10);
maxerr = max(max(abs(V1 - V2)));
data = [];

if opt.Display
   plot(t,V1)
   xlabel('time [\mus]')
   ylabel('V(t)')
   axis tight, grid on, box on
end

end