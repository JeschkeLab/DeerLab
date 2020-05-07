function [pass,maxerr] = test(opt)

% Test manual zero-time correction

t = -5:1:80;
x = 1000 - (t.^2);
x = x/max(x);
t = t + abs(min(t));
ztin = abs(min(t));

[tcorr,ztout] = correctzerotime(x,t,ztin);

% Pass 1: corrected time-axis is equal to original
pass(1) = all(abs(tcorr - t.') < 1e-10);
% Pass 2: zero-time is returned properly
pass(2) = abs(ztout' - ztin) < 1e-10;

pass = all(pass);

maxerr = max(abs(tcorr - t.')); 

end