%
% CHANGEGRAIN Conversion from non-uniform grain to uniform grain
%
%   [Snew,tnew] = CHANGEGRAIN(S,t)
%   Interpolates the input signal (S) on a time axis (t) to a time grain 
%   of 4ns. The resulting signal (Snew) and time axis (tnew) are returned 
%   as arguments. 
%
%   [Snew,tnew] = CHANGEGRAIN(S,t,dtnew)
%   The ouput time grain can be specified by the (dtnew) argument.
%
% Adapted from Gunnar Jeschke by Luis Fabregas.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function [InterpSignal,InterpTimeAxis] = changegrain(Signal,TimeAxis,NewTimeGrain)

max_points = 2048; % maximum number of data points
DefaultTimeGrain = 4; % minimum time increment, should coincide with handles.time_grain in DeerLab.m


OldGrain = TimeAxis(end)/length(TimeAxis);
[IndirectDimension,DirectDimension]=size(Signal);

% If new grain not provided check whether to use the default value
if nargin<2
    NewTimeGrain = OldGrain;
    if NewTimeGrain < DefaultTimeGrain
          NewTimeGrain = DefaultTimeGrain; 
    end
    if mod(NewTimeGrain,DefaultTimeGrain) ~= 0
        NewTimeGrain = DefaultTimeGrain*ceil(NewTimeGrain/DefaultTimeGrain);
    end
end

%Get data reduction factor the change
DataChangeFactor = NewTimeGrain/OldGrain;

max_nexp=max_points; % no more than the specified maximum number of data points

if DirectDimension>max_nexp
    redfac=ceil(DirectDimension/max_nexp);
    if redfac>DataChangeFactor % requested dt is too small
        NewTimeGrain=redfac*OldGrain;
    end
    if mod(NewTimeGrain,DefaultTimeGrain)~=0 % next multiple of eight
        NewTimeGrain=DefaultTimeGrain*ceil(NewTimeGrain/DefaultTimeGrain);
    end
    redfac=round(NewTimeGrain/OldGrain);
    % Data reduction with smoothing
    pas=pascal(redfac);
    b=zeros(1,redfac);
    for k=1:redfac
        b(k)=pas(redfac-k+1,k);
    end
    b=b/sum(b);
    z0=Signal;
    a=1;
    [IndirectDimension,n2]=size(z0);
    Data2Interpolate=zeros(IndirectDimension,n2);
    for kk=1:IndirectDimension
        z0a=z0(kk,:);
        z0a=filter(b,a,z0a);
        Data2Interpolate(kk,:)=z0a;
        Data2Interpolate(kk,1)=z0(kk,1);
    end
else
    Data2Interpolate = Signal;
end

% Make new time axis
StartTime = TimeAxis(1); 
EndTime = TimeAxis(end);
minTime = NewTimeGrain*ceil(StartTime/NewTimeGrain);
if minTime < StartTime
 minTime = StartTime; 
end
maxTime = NewTimeGrain*floor(EndTime/NewTimeGrain);
% while maxt>te, maxt=maxt-dt; end;
if maxTime > EndTime
 maxTime = EndTime; 
end
newLength = 1 + floor((maxTime - minTime)/NewTimeGrain);
maxTime = minTime + NewTimeGrain*(newLength - 1);
InterpTimeAxis = minTime:NewTimeGrain:maxTime;

% Interpolation onto new time axis
InterpSignal = zeros(IndirectDimension,length(InterpTimeAxis));
[IndirectDimension,~]=size(Data2Interpolate);
for kk=1:IndirectDimension
    InterpSignal(kk,:) = interp1(TimeAxis,Data2Interpolate(kk,:),InterpTimeAxis,'pchip',Data2Interpolate(kk,1));
end