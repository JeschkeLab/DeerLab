function [pass,maxerr] = test(opt)

N = 100;
dt = 0.008;
t = linspace(0,dt*N,N);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
S = K*P;

try
    S2 = rand(10,1);
    selregparam(S2,K,r,'tikhonov',{'aic','gml','gcv'});
    err(1) = true;
catch
    err(1) = false;
end

try
    selregparam(S,K,r,'tikhonov','aic','GlobalWeights',[0.2 0.8]);
    err(2) = true;
catch
    err(2) = false;
end


try
    selregparam({S,S},{K,K},r,'tikhonov','aic','GlobalWeights',[1 0.5]);
    err(3) = true;
catch
    err(3) = false;
end

try
    selregparam({S,S},{K,K},r,'tikhonov','aic','NoiseLevel',[0.5]);
    err(3) = true;
catch
    err(3) = false;
end

try
    selregparam({S},{K,K},r,'tikhonov','aic','GlobalWeights',[1 0.5]);
    err(4) = true;
catch
    err(4) = false;
end


try
    selregparam({S,S},{K,K},r,'tikhonov','aic','GlobalWeights',[1 0.5]);
    err(5) = true;
catch
    err(5) = false;
end

try
    S2 = S + 1i*S;
    selregparam({S2,S2},{K,K},r,'tikhonov','aic','GlobalWeights',[1 0.5]);
    err(5) = true;
catch
    err(5) = false;
end

pass = all(err);
 
maxerr = 0;



end