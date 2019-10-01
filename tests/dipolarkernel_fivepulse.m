
function [err,data,maxerr] = test(opt,olddata)

N = 100;
t = linspace(0,3,N);
r = time2dist(t);

K1 = dipolarkernel(t,r,'FivePulseCoeff',[0.5]);
K2 = dipolarkernel(t,r,'FivePulseCoeff',[0.5 max(t)/2]);

err = any(abs(K1 - K2)>1e-10);
maxerr = max(max(abs(K1 - K2)));
data = [];

end