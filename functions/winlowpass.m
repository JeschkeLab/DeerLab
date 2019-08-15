function Signal = winlowpass(Signal,PassBand,SamplingRate,varargin)

[Steepness,MinimalAttenuation] = parseoptional({'Steepness','MinimalAttenuation'},varargin);


if isempty(Steepness)
    Steepness  = 0.85;
else
    validateattributes(Steepness,{'numeric'},{'nonnegative','scalar'},mfilename,'Steepness')
end

if isempty(MinimalAttenuation)
    MinimalAttenuation = 60;
else
    validateattributes(MinimalAttenuation,{'numeric'},{'nonnegative','scalar'},mfilename,'MinimalAttenuation')
end

if PassBand>SamplingRate
   error('Input pass band frequency cannot exceed the sampling rate.') 
end

N = length(Signal);

%Prepare the angular frequencies
PassBand = PassBand*2*pi;
Bandwidth = SamplingRate*2*pi;

%Compute width of transition band as in MATLAB's lowpass function
TransitionBand = (1 - Steepness)*(Bandwidth - PassBand);
Stopband = PassBand + TransitionBand;
TransitionBand = TransitionBand/Bandwidth;

%Normalize all frequencies with respect to bandwith
Stopband = Stopband/Bandwidth*pi;
PassBand  = PassBand/Bandwidth*pi;

%Determine filter order by Kaiser's empirical formula
FilterOrder = ceil((MinimalAttenuation-7.95)/(14.36*TransitionBand/(2*pi))+1) + 1;
FilterAxis = 0:FilterOrder-1;

%Determine Kaiser's beta parameter via his empirical formula
beta = 0.1102*(MinimalAttenuation-8.7);

%Determine cutoff frequency
CutoffBand = (Stopband+PassBand)/2;

%Simulate an ideal IIR low-pass filter
alpha = (FilterOrder-1)/2;
IdealFilterAxis = FilterAxis - alpha + eps;
IIRtransferFcn = sin(CutoffBand*IdealFilterAxis)./(pi*IdealFilterAxis);

%Calculate the Kaiser window function
Window = (kaiser_(FilterOrder,beta))';

%Window the IIR filter to make it FIR
FIRtransferFcn = IIRtransferFcn.*Window;

%Apply FIR low-pass filter to the input signal
[Signal,DelayedSignal] = filter(FIRtransferFcn,1,Signal);

%Calculate FIR filter group delay
Delay = (FilterOrder-1)/2;

%Correct for FIR filter delay
Signal = [Signal DelayedSignal'];
Signal = Signal(Delay+1:Delay+N);

end