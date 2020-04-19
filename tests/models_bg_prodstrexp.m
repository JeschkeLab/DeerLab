function [pass,maxerr] = test(opt)

% Basic functionality test on bg_strexp

model = @bg_prodstrexp;
info = model();

t = linspace(-5,5,500);
par0 = [info.parameters(:).default];
bounds = [info.parameters(:).range];
lower = bounds(1:2:end);
upper = bounds(2:2:end);

B1 = model(t,par0);
B2 = model(t.',par0);
B3 = model(t,lower);
B4 = model(t,upper);
 
t0 = 2.5;
B5 = model(t0,par0);
Bval0 = 2.510999155743982e-08;

% Pass 1: dimensionality is correct
pass(1) = isequal(B1,B2);
pass(2) = iscolumn(B1) && iscolumn(B2);
% Pass 3: there are no NaN values
pass(3) = all(~isnan(B1)) & all(~isnan(B2)) & all(~isnan(B3)) & all(~isnan(B4));
% Pass 2: specific value is reproducible
pass(4) = abs(B5-Bval0) < 1e-8;

pass = all(pass);

maxerr = abs(B5-Bval0);

end
