function [tdipn,trace0,trace1,stddev]=compare_modnorm(tdip0,deer0,tdip1,deer1)
%
% Accepts two DEER time-domain traces deer1 and deer0 that are normalized at t=0 (to V0)
% Scales the modulated part of deer1 so that it has minimum r.m.s.
% deviation from the modulated part of deer0
%
% the time axes of both decays must have the same spacing
% the shorter trace determines the number of points used in r.m.s. error
% computation
%
% Output: unchanged or truncated trace0, scaled and possibly truncated trace1
% and standard deviation of points between traces (stddev)
%
% Used for comparison of DEER decays
%
% (c) G. Jeschke, 2003

dim=size(deer0);
tdipn=tdip0;
if length(deer1)<length(deer0), 
    dim=size(deer1); 
    tdipn=tdip0(1:length(deer1));
end;
ltrc0=log(deer0(1:dim(1),1:dim(2)));
ltrc1=log(deer1(1:dim(1),1:dim(2)));
n=length(ltrc0);

sc=sum(ltrc0.*ltrc0)/sum(ltrc0.*ltrc1);
diff=ltrc0-sc*ltrc1;
stddev=sqrt(sum(diff.*diff)/(n-1));

trace1=exp(sc*ltrc1);
trace0=exp(ltrc0);
