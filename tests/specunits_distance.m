function [err,data,maxerr] = test(opt,olddata)

t = linspace(-0.5,5,100);

r_nm = linspace(0,1.2,100);
r_A = linspace(0,12,100);

%Call some functions which uses units heuristics
K_nm = dipolarkernel(t,r_nm);
K_A = dipolarkernel(t,r_A);

%Should not be equal
err(1) = isequal(K_nm,K_A);

r_nm = specunits(r_nm,'nm');
r_A = specunits(r_A,'A');


%Call some functions which uses units heuristics
K_nm = dipolarkernel(t,r_nm);
K_A = dipolarkernel(t,r_A);

%Should not be equal
err(2) = ~isequal(K_nm,K_A);

err = any(err);
maxerr = 0;
data = [];


end
