function [tdipn,trace0,trace1,stddev]=compare_modnorm(tdip0,deer0,tdip1,deer1)
%
% Accepts two DEER time-domain traces deer1 and deer0 that are normalized at t=0 (to V0)
% Scales the modulated part of deer1 so that it has minimum r.m.s.
% deviation from the modulated part of deer0
% traces must be one-dimensional
%
% the time axis for computation of the scaling factor is the overlapping 
% range of both time axes, the time increment the larger of the two time
% increments
%
% Output: unchanged trace0 (for backward compatibility), scaled trace1 and
% time axis tdipn for which comparison was performed
% and standard deviation of points between traces (stddev)
%
% Used for comparison of DEER decays
%
% (c) G. Jeschke, 2012

% Generate standard time axis
tmin=max([min(tdip0),min(tdip1)]);
tmax=min([max(tdip0),max(tdip1)]);
dt=max([tdip0(2)-tdip0(1),tdip1(2)-tdip1(1)]);
tdipn=tmin:dt:tmax;

ideer0=interp1(tdip0,deer0,tdipn,'pchip','extrap'); % trace0 interpolated tp standard time axis
ideer1=interp1(tdip1,deer1,tdipn,'pchip','extrap'); % trace0 interpolated tp standard time axis

ltrc0=log(ideer0);
ltrc1=log(ideer1);
n=length(ltrc0);

sc=sum(ltrc0.*ltrc0)/sum(ltrc0.*ltrc1);
diff=ltrc0-sc*ltrc1;
stddev=sqrt(sum(diff.*diff)/(n-1));

trace1=exp(sc*log(deer1));
trace0=deer0;