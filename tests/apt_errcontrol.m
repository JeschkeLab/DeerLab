function [pass,maxerr] = test(opt)

% Check error control of apt() towards wrong inputs

Dimension = 200;
t = linspace(-0.6,2,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
S = K*P;

% Pass 1: More than one DDS value passed
aptK = aptkernel(t);
try
    apt(S,aptK,[0.05 125]);
    pass(1) = false;
catch
    pass(1) = true;
end

% Pass 2: APT kernel passed as matrix
K = rand(100,100);
try
    apt(S,K,2);
    pass(2) = false;
catch
    pass(2) = true;
end

% Pass 3: Signal is complex
S = S + 1i*S;
try
    apt(S,aptK,2);
    pass(3) = false;
catch
    pass(3) = true;
end

pass = all(pass);

maxerr = NaN;

end