function [err,data,maxerr] = test(opt,olddata)

r = linspace(3,9,100);

t_us = linspace(-0.5,50,100);
t_ns = linspace(-500,50000,100);

%Call some functions which uses units heuristics
K_us = dipolarkernel(t_us,r);
K_ns = dipolarkernel(t_ns,r);

%Should not be equal
err(1) = isequal(K_us,K_ns);

t_us = specunits(t_us,'us');
t_ns = specunits(t_us,'ns');


%Call some functions which uses units heuristics
K_us = dipolarkernel(t_us,r);
K_ns = dipolarkernel(t_ns,r);

%Should not be equal
err(2) = ~isequal(K_us,K_ns);

err = any(err);
maxerr = 0;
data = [];


end
