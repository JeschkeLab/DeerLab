function [err,data,maxerr] = test(opt,olddata)

%us/nm
t = linspace(0,5,200);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.1]);
S1 = dipolarsignal(t,3);

%us/nm
t = linspace(0,5,200);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.1]);
S2 = dipolarsignal(t,3,'moddepth',1);

S3 = dipolarkernel(t,3);


err = any(abs(S1 - S3)>1e-12);
maxerr = max(max(abs(S1-S3)));
data = [];

if opt.Display
    figure(8),clf
    hold on
    plot(t,S1)
    plot(t,S2)
    plot(t,S3)
end

end