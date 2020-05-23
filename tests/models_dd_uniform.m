function [pass,maxerr] = test(opt)

% Check dimensionality, non-negativity, and values at boundaries of model

info = dd_uniform();

r = linspace(0,30,301);
par0 = [info.parameters(:).default];
bounds = [info.parameters(:).range];
lower = bounds(1:2:end);
upper = bounds(2:2:end);

P1 = dd_uniform(r,par0);
P2 = dd_uniform(r.',par0);
P3 = dd_uniform(r,lower);
P4 = dd_uniform(r,upper);

% Pass 1: dimensionality is correct
pass(1) = isequal(P1,P2);
% Pass 2: non-negativity of default values fulfilled
pass(2) = all(P1 >= 0);
% Pass 3: non-negativity of default boundaries fulfilled
pass(3) = all(P1 >= 0) & all(P2 >= 0);
% Pass 4: there are no NaN values
pass(4) = all(~isnan(P1)) & all(~isnan(P2)) & all(~isnan(P3)) & all(~isnan(P4));
% Pass 5: integral is 1
dr = mean(diff(r));
pass(5) = abs(sum(P1)*dr-1)<1e-10;
% Pass 6: check range
rL = par0(1);
rR = par0(2);
pass(6) = all(P1(r<rL)==0) && all(P1(r>rR)==0);
% Pass 7: compatible with non-uniform distance vectors
rnus = sqrt(linspace(1.5,6^2,800));
P5 = dd_uniform(rnus,par0);
pass(7) = round(trapz(rnus,P5),4) == 1;

pass = all(pass);

maxerr = NaN;

end