%
% CORRECTPHASE Phase correction of complex-valued data

%   [Vr,Vi,ph,io] = CORRECTPHASE(___)
%   ___ = CORRECTPHASE(V)
%   ___ = CORRECTPHASE(V,phase)
%   ___ = CORRECTPHASE(V,phase,true/false)
%   Performs a phase optimization on the complex-valued data (V) by
%   minimization of the imaginary component of the data. The phase can be corrected 
%   manually by specifying a phase (phase), in radians. A third boolean argument
%   can be passed to enable/disable the fitting of a possible offset on the
%   imaginary component of the data. Defaults to false.
%   The function returns the real (Vr) and imaginary (Vi) parts the phase corrected
%   data, the optimized phase (ph) and imaginary offset (io) employed for the correction.
%
%   ___ = CORRECTPHASE(V2D)
%   ___ = CORRECTPHASE(V2D,phases)
%   ___ = CORRECTPHASE(V2D,phases,true/false)
%   Two-dimensional datasets (V2D), e.g. from multiple scans measurements,
%   can be provided, and the phase correction will be done on each trace individually.
%   The first dimension V2D(:,i) must contain the single traces. An array of phases 
%   (phases) can be specified to manually correct the traces.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function [Vreal,Vimag,Phase,ImagOffset] = correctphase(V,Phase,fitImagOffset)


% Input parsing
%-------------------------------------------------------------------------------
switch nargin
    case 0
        error('At least one input (V) is requried.')
    case 1
        Phase = [];
        fitImagOffset = false;
    case 2
        fitImagOffset = false;
    case 3
    otherwise
        error('At most three inputs (V, Phase, fitImagOffset) are possible.');
end

validateattributes(V,{'numeric'},{'2d'},mfilename,'V')
if isvector(V)
    Ntraces = 1;
    V = V(:); % ensure column vector
else
    Ntraces = size(V,2);
end
n = size(V,1);

if ~isempty(Phase)
    validateattributes(Phase,{'numeric'},{'nonnegative'},mfilename,'Phase')
    if numel(Phase)~=Ntraces
       error('The number of input phases must agree with the number of traces.') 
    end
end

validateattributes(fitImagOffset,{'logical'},{'nonempty'},mfilename,'FittedImaginaryOffset')
fitPhase = isempty(Phase);

% Phase/offset fitting
%-------------------------------------------------------------------------------
options = optimset('MaxFunEvals',1e5,'MaxIter',1e5);

ImagOffset = zeros(1,Ntraces);
if fitPhase
    Phase = zeros(1,Ntraces);
    
    for i = 1:Ntraces
        FitRange = round(n/8):n; % use only last 7/8 of data for phase/offset correction
        V_ = V(FitRange,i);
        par0(1) = mean(angle(V_)); % use average phase as initial value
        if fitImagOffset
            par0(2) = mean(imag(V_)); % use average offset as initial value
        end
        fun = @(par)imaginarynorm(par,V_);
        pars = fminsearch(fun,par0,options);
        Phase(i) = pars(1);
        if fitImagOffset
            ImagOffset(i) = pars(2);
        end
    end

else
    if fitImagOffset
        % Fit only imaginary offset        
        for i = 1:Ntraces
            par0 = 0;
            fun = @(offset)imaginarynorm([Phase offset],V(:,i));
            ImagOffset(i) = fminsearch(fun,par0);
        end
    end
end


ImagOffset = ImagOffset*1i;

% Apply phase/offset correction
ph = exp(1i*Phase);
Vc = (V - ImagOffset)./ph;

% Output
Vreal = real(Vc);
Vimag = imag(Vc);
Phase = angle(ph); % map phase angle to [-pi,pi) interval

end

function ImagNorm = imaginarynorm(params,V)
% Computes norm of the imaginary part of phase-corrected data from zero before
% phase correction, an offset can be subtracted from the imaginary part.

phase = params(1);
if numel(params)>1
    imoffsets = params(2);
else
    imoffsets = 0;
end

Vcorr = (V-1i*imoffsets).*exp(-1i*phase);
ImagNorm = norm(imag(Vcorr));

end
