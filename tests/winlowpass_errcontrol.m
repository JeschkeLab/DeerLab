function [pass,maxerr] = test(opt)

dt = 0.5e-9;
nu1 = 0.5e9;
N = 300;
t = linspace(0,dt*N,N);
S = exp(-5e7*t).*(cos(2*pi*nu1*t));
sampl = 1/2*(1/dt);
wp = 0.01e9;
ws =  0.04e9;

try
    winlowpass(S,wp,ws,sampl);
    err(1) = true;
catch
    err(1) = false;
end

try
    winlowpass(S,ws,sampl,wp);
    err(2) = true;
catch
    err(2) = false;
end


pass = all(err);
 
maxerr = 0;


end