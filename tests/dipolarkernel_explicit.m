function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

load(fullfile('comparison','oldDAkernel500'));

N = 500;
dt = 0.008;
t = linspace(0,dt*N,N);
DistAxis = time2dist(t);
kernelOut = dipolarkernel(t,DistAxis,'KCalcMethod','explicit');
kernelOut = kernelOut/mean(diff(DistAxis));
err = any(abs(kernelOut - kernel)>7e-3);
maxerr = max(max(abs(kernelOut - kernel)));
data = [];


if opt.Display
    figure(8),clf
    subplot(121)
    imagesc(kernelOut)
    subplot(122)
    hold on
%     plot(t,kernel(:,5))
    plot(t,kernelOut(:,5))
end

end