function [pass,maxerr] = test(opt)

% Check indifference of winlowpass() towards input dimensionality

t = linspace(0,5,80);
S = dipolarsignal(t,3);

Sfilt1 = winlowpass(S,3e6,2e6,8e6);
Sfilt2 = winlowpass(S.',3e6,2e6,8e6);

% Pass 1:  both noise vectors are column vectors
pass(1) = iscolumn(Sfilt2) & iscolumn(Sfilt1);
% Pass 2: both noise vectors are equal
pass(2) = isequal(Sfilt2,Sfilt1);

pass = all(pass);

maxerr = NaN; 

end