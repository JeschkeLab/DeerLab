function [tdipn,trace0,trace1,stddev]=compare_modnorm_linear(tdip0,deer0,tdip1,deer1,depth0,depth1)
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
% Used for comparison of form factors
%
% (c) G. Jeschke, 2012

trace0=deer0;

% Generate standard time axis
tmin=max([min(tdip0),min(tdip1)]);
tmax=min([max(tdip0),max(tdip1)]);
dt=max([tdip0(2)-tdip0(1),tdip1(2)-tdip1(1)]);
tdipn=tmin:dt:tmax;

ideer0=interp1(tdip0,deer0,tdipn,'pchip','extrap')-1+depth0; % trace0 interpolated tp standard time axis
ideer1=interp1(tdip1,deer1,tdipn,'pchip','extrap')-1+depth1; % trace0 interpolated tp standard time axis

n=length(ideer0);

ideer1p=depth0*(ideer1-1+depth1)/depth1+1-depth0;
diff=ideer0-ideer1p;
stddev=sqrt(sum(diff.*diff)/(n-1));


trace1=depth0*(deer1-1+depth1)/depth1+1-depth0;
