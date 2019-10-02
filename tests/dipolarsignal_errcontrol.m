function [err,data,maxerr] = test(opt,olddata)

N = 100;
t = linspace(0,3,N);
r = time2dist(t);
P = rd_onegaussian(r,[4,0.4]);

try
    dipolarsignal(t,r);
    err(1) = true;
catch
    err(1) = false;
end

try
    dipolarsignal(t,r,P,'moddepth',-0.1);
    err(2) = true;
catch
    err(2) = false;
end

try
    dipolarsignal(t,r,P,'moddepth',1.5);
    err(3) = true;
catch
    err(3) = false;
end

try
    r2 = r.^2;
    dipolarsignal(t,r2,P);
    err(4) = true;
catch
    err(4) = false;
end

try
    r2 = linspace(2,5,10);
    dipolarsignal(t,r2,P);
    err(5) = true;
catch
    err(5) = false;
end

err = any(err);
data = [];
maxerr = 0;

end