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

tstart1 =  backgroundstart(S,t,@td_poly1);
tstart2 =  backgroundstart(S,t.',@td_poly1);

tstart3 =  backgroundstart(S,t,@td_poly1);
tstart4 =  backgroundstart(S.',t,@td_poly1);


err(1) = abs(tstart1 - tstart2)>1e-10;
err(2) = abs(tstart3 - tstart4)>1e-10;
err = any(err);
data = [];
maxerr = abs(tstart1 - tstart2);


end