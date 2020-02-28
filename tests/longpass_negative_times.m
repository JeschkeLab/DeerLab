function [pass,maxerr] = test(opt)

Dimension = 100;
t = linspace(-0.5,5,Dimension);
r = time2dist(t);
P = rd_twogaussian(r,[4,0.3,6.5,0.3,0.5]);
K = dipolarkernel(t,r);
DipEvoFcn = K*P;

Filtered = longpass(t,DipEvoFcn,2);

Result = fitregmodel(Filtered,K,r,'tikhonov','aic','Solver','fnnls');

error = abs(Result - P);
err(1) = any(error>7e-1);
maxerr = max(error);
pass = all(err);
 

if opt.Display
    figure(9),clf
    subplot(132),hold on
    plot(r,P)
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