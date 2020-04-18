function [pass,maxerr] = test(opt)

% Basic functionality test on bg_exvol

% Pass 1: dimensionality is correct
t = linspace(0,6,301);
R = 14; % nm
lam = 1;
c = 500; % uM
par = [R lam*c];
B1 = bg_exvol(t,par);
B2 = bg_exvol(t.',par);
pass(1) = iscolumn(B1) && iscolumn(B2);

% Pass 2: specific value
t = 6; % us
R = 14; % nm
lam = 1;
c = 500; % uM
Bval = bg_exvol(t,[R lam*c]);
Bval0 = 0.5011333;
pass(2) = abs(Bval-Bval0);

pass = all(pass);

maxerr = abs(Bval-Bval0);

end
