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
%   ct = FFTSPEC(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order. 
%
%   'Type' - Type of spectrum to be returned ('real','imag','abs')
%
%   'ZeroFilling' - Number of elements in the output FFT spectrum 
%                   (default = 2*length(S)).
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

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