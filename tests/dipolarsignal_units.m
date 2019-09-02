function [err,data,maxerr] = test(opt,olddata)

%us/nm
t = linspace(0,5,200);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.1]);
S1 = dipolarsignal(t,r,P);

%us/A
t = linspace(0,5,200);
r = 10*time2dist(t);
S2 = dipolarsignal(t,r,P);

%ns/A
t = 1000*linspace(0,5,200);
r = 10*time2dist(t);
S3 = dipolarsignal(t,r,P);

%ns/nm
t = 1000*linspace(0,5,200);
r = time2dist(t);
S4 = dipolarsignal(t,r,P);

err = isequal(S1,S2,S3,S4);
maxerr = max(max(abs(S1-S2)));
data = [];

if opt.Display
    figure(8),clf
    hold on
    plot(t,S1)
    plot(t,S2)
    plot(t,S3)
    plot(t,S4)
    legend('us/nm','us/A','ns/A','ns/nm')
end

end