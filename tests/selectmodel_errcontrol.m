function [err,data,maxerr] = test(opt,olddata)

N = 200;
dt = 0.016;
t = linspace(0,dt*N,N);
r = time2dist(t);
InputParam = [3 0.3 5 0.3 0.5];
P = rd_twogaussian(r,InputParam);
K = dipolarkernel(t,r);
S = K*P;
Models = {@rd_onegaussian,@rd_twogaussian,@rd_threegaussian};

try
    selectmodel(Models,S,r,K,'gsad'); 
    err(1) = true;
catch
    err(1) = false;
end
try
    selectmodel(Models); 
    err(2) = true;
catch
    err(2) = false;
end
try
    selectmodel(Models,S,r,K,'gsad'); 
    err(3) = true;
catch
    err(3) = false;
end

try
    selectmodel(Models,S,r,K,'aic','upper',[1 2 3]); 
    err(3) = true;
catch
    err(3) = false;
end
try
    selectmodel(Models,S,r,K,'aic','lower',[1 2 3]); 
    err(4) = true;
catch
    err(4) = false;
end
try
    selectmodel(Models,K,'aic','lower',{[1 2 3],[1 2 3],[1 3 4],[1 4 5]}); 
    err(5) = true;
catch
    err(5) = false;
end
try
    selectmodel(@rd_fivegaussian,K,'aic','lower',{[1 2 3],[1 2 3],[1 3 4],[1 4 5]}); 
    err(6) = true;
catch
    err(6) = false;
end
try
    selectmodel(Models,S,r,K,'aic',[1 2 3]); 
    err(7) = true;
catch
    err(7) = false;
end

err = any(err);
maxerr = NaN;
data = [];

end

