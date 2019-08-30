function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

Ntime = 302;
Ndist = 200;

TimeStep = 0.008;
TimeAxis = linspace(-TimeStep*Ntime/2,TimeStep*Ntime/2,Ntime);
[~,rmin,rmax] = time2dist(TimeAxis);
DistAxis = linspace(rmin,rmax,Ndist);
kernel = dipolarkernel(TimeAxis,DistAxis,'KernelCalcMethod','explicit');
kernel = kernel/mean(diff(DistAxis));
trace = kernel(:,1);

negtrace = trace(1:Ntime/2)';
postrace = trace(Ntime/2+1:end)';

err(1) = any(any(abs(fliplr(negtrace) - postrace)>1e-12));
err = any(err);
maxerr = max(abs(fliplr(negtrace) - postrace));
data = [];

if opt.Display
    figure(8),clf
    subplot(121)
    imagesc(kernel)
    subplot(122)
    hold on
    plot(TimeAxis(Ntime/2+1:end),fliplr(negtrace))
    plot(TimeAxis(Ntime/2+1:end),postrace)
end
end