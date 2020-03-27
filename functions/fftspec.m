%
% FFTSPEC Fast-Fourier transform spectrum
%
%   spec = FFTSPEC(t,S)
%   Computes the magnitude FFT spectrum of the signal (S) on the time axis
%   (t).
%
%   [nu,spec] = FFTSPEC(t,S)
%   If two output arguments are requested, the frequency axis is returned
%   as well.
%
%   [nu,spec] = FFTSPEC(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order.
%
%   'Type' - Type of spectrum to be returned ('real','imag','abs')
%
%   'ZeroFilling' - Number of elements in the output FFT spectrum
%                   (default = 2*length(S)).
%
%   'Apodization' - Use apodization window
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function varargout = fftspec(t,S,varargin)

if nargin<2
    error('Not enough inputs.');
end

%Parse optional input
[Type,ZeroFilling,useApodization] = parseoptional({'Type','ZeroFilling','Apodization'},varargin);

if isempty(Type)
    Type = 'abs';
else
    validateattributes(Type,{'char'},{},mfilename,'Type')
    allowedInput = {'abs','real','imag','complex'};
    validatestring(Type,allowedInput);
end

if isempty(useApodization)
    useApodization = true;
end
if isempty(ZeroFilling)
    ZeroFilling = 2*length(S);
else
    validateattributes(ZeroFilling,{'numeric'},{'scalar','nonnegative'},mfilename,'ZeroFilling')
end
validateattributes(useApodization,{'logical'},{'nonempty'},mfilename,'useApodization')
validateattributes(S,{'numeric'},{'2d','nonempty'},mfilename,'FitData')
validateattributes(t,{'numeric'},{'2d','nonempty','increasing'},mfilename,'t')
%Use column vectors
t = t(:);
S = S(:);

%If requested apply Hamming apodization window
if useApodization
    arg = linspace(0,pi,length(S));
    arg = arg.';
    ApoWindow = 0.54 + 0.46*cos(arg);
    S = S.*ApoWindow;
end

%Compute fft spectrum
Spectrum = fftshift(fft(S,ZeroFilling));

%Get the requested component/type of spectrum
switch Type
    case 'abs'
        Spectrum = abs(Spectrum);
    case 'real'
        Spectrum  = real(Spectrum);
    case 'imag'
        Spectrum = imag(Spectrum);
end

%If requested, get frequency axis and switch output order
if nargout>1
    FreqAxis = time2freq(t,ZeroFilling);
    varargout{1} = FreqAxis;
    varargout{2} = Spectrum;
else
    varargout{1} = Spectrum;
end


end