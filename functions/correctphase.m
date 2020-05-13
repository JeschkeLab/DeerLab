%
% CORRECTPHASE Phase correction of complex-valued data

%   [Vr,Vi,ph,io] = CORRECTPHASE(___)
%   ___ = CORRECTPHASE(V)
%   ___ = CORRECTPHASE(V,phase)
%   ___ = CORRECTPHASE(V,phase,true/false)
%   Performs a phase optimization on the complex-valued data (V) by
%   minimization of the imaginary component of the data. The phase can be corrected 
%   manually by specifying a phase (phase), in radians. A third boolean argument
%   can be passed to enable/diasable the fitting of a possible offset on the
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

%--------------------------------------------------------------------------
%Input parsing
%--------------------------------------------------------------------------
switch nargin
    case 1
        Phase = [];
        fitImagOffset = false;
    case 2
        fitImagOffset = false;
    case 3
    otherwise
        error('Wrong number of input arguments.');
end

if all(size(V)>1)
    Ntraces = size(V,2);
else
    Ntraces = 1;
    %Ensure column vector
    V = V(:);
end

validateattributes(V,{'numeric'},{},mfilename,'PrimaryData')
if ~isempty(Phase)
    validateattributes(Phase,{'numeric'},{'nonnegative'},mfilename,'Phase')
    if numel(Phase)~=Ntraces
       error('The number of input phases must agree with the number of traces.') 
    end
end
validateattributes(fitImagOffset,{'logical'},{'nonempty'},mfilename,'FittedImaginaryOffset')

options = optimset('MaxFunEvals',1e5,'MaxIter',1e5);

ImagOffset = zeros(1,Ntraces);
% If phase is not provided, then fit it (and imag. offset)
if isempty(Phase)
    for i=1:Ntraces
        phi0 = mean(angle(V(:,i))); % use average phase as starting point for fit
        pars(1) = phi0;
        if fitImagOffset
            pars(2) = 0;
        end
        FitStart = round(size(V,1)/8); % use only last 7/8 of data for phase/offset correction
        pars = fminsearch(@(par)imaginarynorm(par,V(FitStart:end,i)),pars,options);
        Phase(i) = pars(1);
        if fitImagOffset
            ImagOffset(i) = pars(2);
        end
    end

elseif ~isempty(Phase) && fitImagOffset
    for i=1:Ntraces
        ImagOffset(i) = fminsearch(@(offset)imaginarynorm([Phase offset],V(:,i)),0);
    end
end

%Wrap phases to [0 pi] range
Phase = mod(Phase,pi);

%Do phase correction and normalize
Vc = (V - 1i*ImagOffset).*exp(1i*Phase);
Vreal = real(Vc);
Vimag = imag(Vc);

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
if any(imoffsets<0)
    ImagNorm = 1e6;
    return
end

ImagComponent = imag((V-1i*imoffsets).*exp(1i*phase));
ImagNorm = norm(ImagComponent);

end