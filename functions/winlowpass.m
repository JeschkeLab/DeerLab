%
% WINLOWPASS Windowed low-pass FIR filter
%
%   X = WINLOWPASS(S,Fstop,Fpass,Fsamp)
%   Filters the N-point signal S using a Kaiser-windowed lowpass FIR filter
%   with a transition band given by the passband Fpass and stopband  Fstop
%   frequencies in Hz. The filter is then constructed according to the sampling
%   rate Fsamp of the signal in Hz.
%
%   [X,H] = WINLOWPASS(S,Fstop,Fpass,Fsamp)
%   The filter transfer function parameters H requested as a second ouput
%   argument.
%
%   [X,H] = WINLOWPASS(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order. 
%
%   'MinimalAttenuation' - Minimal attenuation level [dB] of the first 
%                           sidelobe after the stopband
%
%   'ForwardBackward' - Enable/disable forward-backward filtering of the 
%                       signal for zero-phase filtering                    
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [FilteredSignal,FIRtransferFcn] = winlowpass(Signal,StopBand,PassBand,SamplingRate,varargin)

%Parse optional input arguments
[MinimalAttenuation,ForwardBackward] = parseoptional({'MinimalAttenuation','ForwardBackward'},varargin);

%Validate inputs
if isempty(MinimalAttenuation)
    MinimalAttenuation = 50;
else
    validateattributes(MinimalAttenuation,{'numeric'},{'nonnegative','scalar'},mfilename,'MinimalAttenuation')
end

if isempty(ForwardBackward)
    ForwardBackward = true;
else
    validateattributes(ForwardBackward,{'logical'},{'scalar','nonempty'},mfilename,'ForwardBackward')
end

validateattributes(Signal,{'numeric'},{'2d','nonempty'},mfilename,'Signal')
validateattributes(StopBand,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'Signal')
validateattributes(PassBand,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'Signal')
validateattributes(SamplingRate,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'Signal')

if PassBand>SamplingRate
    error('Input pass band frequency cannot exceed the sampling rate.')
end

if StopBand<PassBand
    error('Stopband frequency must be larger than the passband frequency')
end

N = length(Signal);

%Normalize frequencies
StopBand = StopBand/SamplingRate;
PassBand = PassBand/SamplingRate;

%Determine filter order by Kaiser's empirical formula
TransitionBand = StopBand - PassBand;
FilterOrder = ceil((MinimalAttenuation-7.95)/(14.36*TransitionBand)+1) + 1;
FilterAxis = 0:FilterOrder-1;
%Determine Kaiser's beta parameter via his empirical formula
beta = 0.1102*(MinimalAttenuation-8.7);

%Determine cutoff frequency
CutoffBand = (StopBand+PassBand)/2;
%Simulate an ideal IIR low-pass filter
alpha = (FilterOrder-1)/2;
IdealFilterAxis = FilterAxis - alpha + eps;
IIRtransferFcn = sin(2*pi*CutoffBand*IdealFilterAxis)./(pi*IdealFilterAxis);

%Calculate the Kaiser window function
Window = (kaiser_(FilterOrder,beta))';

%Window the IIR filter to make it FIR
FIRtransferFcn = IIRtransferFcn.*Window;

%Normalize the transfer function coefficients
FIRtransferFcn = FIRtransferFcn/sum(FIRtransferFcn);

if ForwardBackward
    %If requested create a forward/backward FIR filter transfer function
    FIRtransferFcn = FIRtransferFcn.*fliplr(FIRtransferFcn);
    %Normalize again
    FIRtransferFcn = FIRtransferFcn/sum(FIRtransferFcn);
end

%Apply FIR low-pass filter to the input signal
FilteredSignal = filter(FIRtransferFcn,1,[Signal zeros(length(FIRtransferFcn),1)'] );

%Calculate FIR filter group delay
Delay = FilterOrder/2 - 1;

%Correct for FIR filter delay
FilteredSignal = FilteredSignal(Delay+1:Delay+N);

end