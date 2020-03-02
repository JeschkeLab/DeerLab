function [pass,maxerr] = test(opt)

% Check error control of dipolarsignal() towards wrong inputs

N = 100;
t = linspace(0,3,N);
r = time2dist(t);
P = rd_onegaussian(r,[4,0.4]);

% Pass 1: not enough input arguments
try
    dipolarsignal(t,r);
    pass(1) = false;
catch
    pass(1) = true;
end

% Pass 2: negative modulation depth
try
    dipolarsignal(t,r,P,'moddepth',-0.1);
    pass(2) = false;
catch
    pass(2) = true;
end

% Pass 3: modulation depth larger than 1
try
    dipolarsignal(t,r,P,'moddepth',1.5);
    pass(3) = false;
catch
    pass(3) = true;
end

pass = all(pass);
 
maxerr = NaN;

end