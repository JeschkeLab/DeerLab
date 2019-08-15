function [err,data,maxerr] = test(opt,olddata)

%======================================================
% Check kernel is constructed properly
%======================================================

load(fullfile('comparison','oldDAkernel500'));

N = 500;
TimeStep = 0.008;
TimeAxis = linspace(0,TimeStep*N,N);
DistAxis = time2dist(TimeAxis);
kernelOut = dipolarkernel(TimeAxis,DistAxis,[],'KernelCalcMethod','explicit');
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
%     plot(TimeAxis,kernel(:,5))
    plot(TimeAxis,kernelOut(:,5))
end

end