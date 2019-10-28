function [err,data,maxerr] = test(opt,olddata)

N = 100;
t = linspace(0,3,N);
r = time2dist(t);

try
K = dipolarkernel(t,r,'interference',[0.5 max(t)/2 1]);
err = true;
catch 
err = false;    
end

r2 = r.^2;
try
K = dipolarkernel(t,r2,'interference',[0.5 max(t)/2 1]);
err = true;
catch 
err = false;    
end


maxerr = 0;
data = [];

end