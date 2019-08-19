function [err,data,maxerr] = test(opt,olddata)

Dimension = 100;
TimeAxis = linspace(-0.5,5,Dimension);
DistanceAxis = time2dist(TimeAxis);
Distribution = 0.5*gaussian(DistanceAxis,4,0.3) + 0.5*gaussian(DistanceAxis,6.5,0.3);
Distribution = Distribution/sum(Distribution)/mean(diff(DistanceAxis));
Kernel = dipolarkernel(TimeAxis,DistanceAxis);
DipEvoFcn = Kernel*Distribution;

Filtered = longpass(TimeAxis,DipEvoFcn,2);

L = regoperator(Dimension,2);
RegParam = regparamrange(Kernel,L);
RegParam2 = selregparam(RegParam,Filtered,Kernel,L,'tikhonov','aic');
Result = fitregmodel(Filtered,Kernel,DistanceAxis,L,'tikhonov',RegParam2,'Solver','fnnls');

error = abs(Result - Distribution);
err(1) = any(error>3e-1);
maxerr = max(error);
err = any(err);
data = [];

if opt.Display
    figure(9),clf
    subplot(132),hold on
    plot(DistanceAxis,Distribution)
    plot(DistanceAxis,Result)
    subplot(131),hold on
    plot(TimeAxis,DipEvoFcn)
    plot(TimeAxis,Filtered)
    subplot(133),hold on
    plot(abs(fftshift(fft(DipEvoFcn(find(TimeAxis==0):end)))))
    plot(abs(fftshift(fft(Filtered(find(TimeAxis==0):end)))))   
    axis tight
end

end