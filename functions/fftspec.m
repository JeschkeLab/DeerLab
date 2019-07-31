function varargout = fftspec(TimeAxis,Signal,varargin)

if nargin<2
    error('Not enough inputs.');
end
[Type,ZeroFilling] = parseoptional({'Type','ZeroFilling'},varargin);
if isempty(Type)
    Type = 'abs';
else
    validateattributes(Type,{'char'},{},mfilename,'Type')
    allowedInput = {'abs','real','imag','complex'};
    validatestring(Type,allowedInput);
end
if isempty(ZeroFilling)
    ZeroFilling = 2*length(Signal);
else
    validateattributes(ZeroFilling,{'numeric'},{'scalar','nonnegative'},mfilename,'ZeroFilling')
end
validateattributes(Signal,{'numeric'},{'2d','nonempty'},mfilename,'FitData')
validateattributes(TimeAxis,{'numeric'},{'2d','nonempty','nonnegative','increasing'},mfilename,'TimeAxis')

Spectrum = fftshift(fft(Signal,ZeroFilling));

switch Type
    case 'abs'
        Spectrum = abs(Spectrum);
    case 'real'
        Spectrum  = real(Spectrum);
    case 'imag'
        Spectrum = imag(Spectrum);    
end


if nargout>1
    FreqAxis = time2freq(TimeAxis,ZeroFilling);
    varargout{1} = FreqAxis;
    varargout{2} = Spectrum;
else
    varargout{1} = Spectrum;
end


end