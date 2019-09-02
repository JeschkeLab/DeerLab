function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

t = linspace(0,5,200);
%nm
r = time2dist(t);
Knm = dipolarkernel(t,r);

%A
r = 10*time2dist(t);
KA = dipolarkernel(t,r);

err = any(any(abs(KA-Knm)>1e-12));
maxerr = max(max(abs(KA-Knm)));
data = [];

if opt.Display
    figure(8),clf
    hold on
    plot(t,Knm(:,1))
    plot(t,KA(:,1))
end
end