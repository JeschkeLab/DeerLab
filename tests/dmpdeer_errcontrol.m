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

try
    prob = [0.8 0.8 0.8];
    td_dmpdeer(t,r,P,taus,ts,prob);
    err(1) = true;
catch
    err(1) = false;
end

try
    taus = [123 123 123];
    ts = [123 123];
    td_dmpdeer(t,r,P,taus,ts,prob);
    err(2) = true;
catch
    err(2) = false;
end

try
    taus = 'wrong';
    td_dmpdeer(t,r,P,taus,ts,prob);
    err(3) = true;
catch
    err(3) = false;
end


err = any(err);
maxerr = NaN;
data = [];

if opt.Display
   plot(t,V1)
   xlabel('time [\mus]')
   ylabel('V(t)')
   axis tight, grid on, box on
end

end