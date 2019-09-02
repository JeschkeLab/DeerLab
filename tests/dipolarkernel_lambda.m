function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

t = linspace(0,5,200);
%nm
r = time2dist(t);
lambda = 0.4;
K = dipolarkernel(t,r,lambda,'gValue',2);

err = false;
maxerr = 0;
data = [];

if opt.Display
    figure(8),clf
    plot(t,K(:,10))
end
end