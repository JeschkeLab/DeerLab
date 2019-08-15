function FilteredSignal = winlowpass(Signal,StopBand,PassBand,SamplingRate,varargin)

[MinimalAttenuation] = parseoptional({'MinimalAttenuation'},varargin);

if isempty(MinimalAttenuation)
    MinimalAttenuation = 50;
else
    validateattributes(MinimalAttenuation,{'numeric'},{'nonnegative','scalar'},mfilename,'MinimalAttenuation')
end

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

%Apply FIR low-pass filter to the input signal
FilteredSignal = filter(FIRtransferFcn,1,[Signal zeros(length(FIRtransferFcn),1)'] );

%Calculate FIR filter group delay
Delay = (FilterOrder-1)/2;

%Correct for FIR filter delay
FilteredSignal = FilteredSignal(Delay+1:Delay+N);

end