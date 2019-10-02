function [err,data,maxerr] = test(opt,olddata)

%==============================================================
% Get start of background fit ensure that integer is returned
%==============================================================
%Parameters
k = 0.5;
N = 200;
dt = 0.016;
t = linspace(0,N*dt,N);
%Construct some dipolar evolution function 
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);
%Construct background
B = exp(-k*t).';
lam0 = 0.5;
%Account modulation depth for the offset=1
S = dipolarsignal(t,r,P,'background',B,'moddepth',lam0);

try
    backgroundstart(S,t);
    err(1) = true;
catch
    err(1) = false;
end


S2 = dipolarsignal(t,r,P,'background',B,'moddepth',lam0,'phase',pi/2);
try
    backgroundstart(S2,t,@td_poly1);
    err(2) = true;
catch
    err(2) = false;
end

model = rand(2,2);
try
    backgroundstart(S,t,model);
    err(3) = true;
catch
    err(3) = false;
end

try
    backgroundstart(S,t,@td_poly1,'RelSearchStart',0.9,'RelSearchEnd',0.1);
    err(4) = true;
catch
    err(4) = false;
end
err = any(err);
data = [];
maxerr = 0;

end