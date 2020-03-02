function [pass,maxerr] = test(opt)

% Check error control of winlowpass() towards wrong inputs

dt = 0.5e-9;
nu1 = 0.5e9;
N = 300;
t = linspace(0,dt*N,N);
S = exp(-5e7*t).*(cos(2*pi*nu1*t));
sampl = 1/2*(1/dt);
wp = 0.01e9;
ws =  0.04e9;

% Pass 1: stopband larger than the passband
try
    winlowpass(S,wp,ws,sampl);
    pass(1) = false;
catch
    pass(1) = true;
end

% Pass 2: sampling frequency too small
try
    winlowpass(S,ws,sampl,wp);
    pass(2) = false;
catch
    pass(2) = true;
end


pass = all(pass);
 
maxerr = 0;


end