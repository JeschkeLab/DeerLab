function [err,data,maxerr] = test(opt,olddata)

clear aptkernel
N = 200;
dt = 0.008;
t = linspace(0,dt*(N-1),N);

try
    K = aptkernel(t,'ExcitationBandwidth','a');
    err(1) = true;
catch
    err(1) = false;
end

err = any(err);
maxerr = 0;
data = [];

end