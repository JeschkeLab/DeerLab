function [pass,maxerr] = test(opt)

% Chech that dipolar kernel with limited excitation bandwidth works well

t = linspace(0,3.2,80);
r = time2dist(t);
P = dd_onegauss(r,[3,0.5]);
%Use a very limited bandwidth
ExcitationBandwidth = 0.05; %MHz
K = dipolarkernel(t,r,'ExcitationBandwidth',ExcitationBandwidth);
S  = round(K*P,10);

% Pass: all signal elements are zero
pass = all(S==0);

maxerr = NaN;
 
end