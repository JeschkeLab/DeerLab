
function [err,data,maxerr] = test(opt,olddata)

N = 100;
t = linspace(0,3,N);
r = time2dist(t);

lam = 0.5;
B = ones(1,numel(t));

K1 = dipolarkernel(t,r,lam,'interference',[0.5 max(t)/4 0.2 max(t)/2]);
K2 = dipolarkernel(t,r,lam,B,'interference',{0.5 max(t)/4 0.2 max(t)/2});

err = any(abs(K1 - K2./B)>1e-10);
maxerr = max(max(abs(K1 - K2)));
data = [];

end