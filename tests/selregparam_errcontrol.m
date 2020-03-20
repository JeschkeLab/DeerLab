function [pass,maxerr] = test(opt)

% Check error control of selregparam() towards wrong inputs

N = 100;
dt = 0.008;
t = linspace(0,dt*N,N);
r = time2dist(t);
P = dd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
S = K*P;

% Pass 1: signal not right size
try
    S2 = rand(10,1);
    selregparam(S2,K,r,'tikhonov',{'aic','gml','gcv'});
    pass(1) = false;
catch
    pass(1) = true;
end

% Pass 2: too many global weights
try
    selregparam(S,K,r,'tikhonov','aic','GlobalWeights',[0.2 0.8]);
    pass(2) = false;
catch
    pass(2) = true;
end


% Pass 3: not enough noise levels
try
    selregparam({S,S},{K,K},r,'tikhonov','aic','NoiseLevel',[0.5]);
    pass(3) = false;
catch
    pass(3) = true;
end

% Pass 4: not enough signals
try
    selregparam({S},{K,K},r,'tikhonov','aic','GlobalWeights',[1 0.5]);
    pass(4) = false;
catch
    pass(4) = true;
end

% Pass 5: a complex-valued signal is passed
try
    S2 = S + 1i*S;
    selregparam({S2,S2},{K,K},r,'tikhonov','aic','GlobalWeights',[1 0.5]);
    pass(5) = false;
catch
    pass(5) = true;
end

pass = all(pass);
 
maxerr = NaN;



end