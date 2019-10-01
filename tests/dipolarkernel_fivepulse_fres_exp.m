
function [err,data,maxerr] = test(opt,olddata)

N = 100;
t = linspace(0,3,N);
r = time2dist(t);

kernelF = dipolarkernel(t,r,'Method','fresnel','FivePulseCoeff',[0.5 max(t)/2]);
kernelE = dipolarkernel(t,r,'Method','explicit','FivePulseCoeff',[0.5 max(t)/2]);
kernelF = kernelF/mean(diff(r));
kernelE = kernelE/mean(diff(r));

err = any(abs(kernelF - kernelE)>1e-3);
maxerr = max(max(abs(kernelF - kernelE)));
data = [];

end