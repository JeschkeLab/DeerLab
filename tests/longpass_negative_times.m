function [err,data,maxerr] = test(opt,olddata)

Dimension = 100;
t = linspace(-0.5,5,Dimension);
r = time2dist(t);
Distribution = rd_twogaussian(r,[4,0.3,6.5,0.3,0.5]);
K = dipolarkernel(t,r);
DipEvoFcn = K*Distribution;

Filtered = longpass(t,DipEvoFcn,2);

L = regoperator(Dimension,2);
RegParam = regparamrange(K,L);
RegParam2 = selregparam(RegParam,Filtered,K,L,'tikhonov','aic');
Result = fitregmodel(Filtered,K,r,L,'tikhonov',RegParam2,'Solver','fnnls');

error = abs(Result - Distribution);
err(1) = any(error>3e-1);
maxerr = max(error);
err = any(err);
data = [];

if opt.Display
    figure(9),clf
    subplot(132),hold on
    plot(r,Distribution)
    plot(r,Result)
    subplot(131),hold on
    plot(t,DipEvoFcn)
    plot(t,Filtered)
    subplot(133),hold on
    plot(abs(fftshift(fft(DipEvoFcn(find(t==0):end)))))
    plot(abs(fftshift(fft(Filtered(find(t==0):end)))))   
    axis tight
end

end