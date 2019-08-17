%
% CORRECTPHASE Phase correction of complex-valued data
%
%   S = CORRECTPHASE(C)
%   Performs a phase optimization on the complex-valued data (C) by
%   minimization of the imaginary component of the data. The phase
%   corrected data (S) is returned normalized.
%
%   S = CORRECTPHASE(C,p)
%   A phase can be passed manually for the correction by passing a second
%   argument (p).
%
%   S = CORRECTPHASE(C,p,true/false)
%   A third boolean argument can be passed to enable/diasable the fitting 
%   of a possible offset on the imaginary component of the data. Defaults
%   to false.
%
%   [S,p,io] = CORRECTPHASE(C)
%   Additional output arguments can be requested to return the optimal
%   phase (p) and imaginary offset (io) employed for the correction.
%
%
% Adapted from Gunnar Jeschke, DeerAnalysis 2018
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function [correctedSignal,Phase,ImagOffset] = correctphase(PrimaryData,Phase,FittedImaginaryOffset)

%--------------------------------------------------------------------------
%Input parsing
%--------------------------------------------------------------------------
if nargin<3
    FittedImaginaryOffset = false;
else
    validateattributes(FittedImaginaryOffset,{'logical'},{'nonempty'},mfilename,'FittedImaginaryOffset')
end

if iscolumn(PrimaryData)
   PrimaryData = PrimaryData'; 
end
validateattributes(PrimaryData,{'numeric'},{},mfilename,'PrimaryData')

if nargin==1 || nargin==3
    Phase = [];
else
    validateattributes(Phase,{'numeric'},{'scalar','nonnegative'},mfilename,'Phase')
end
DirectDimension = size(PrimaryData,1);

if DirectDimension>1
    Signal=PrimaryData(1,:);
else
    Signal=PrimaryData;
end

% If phase is not provided, then fit it
if isempty(Phase)
    SignalEnd = PrimaryData(length(PrimaryData));
    phi0 = atan2(imag(SignalEnd),real(SignalEnd));
    FittedPhase = phi0;
    FitStart = round(length(Signal)/8); % use only last 7/8 of data for phase/offset correction
    if ~FittedImaginaryOffset
        FittedPhase = FittedPhase(1);
    else
        FittedPhase(2) = 0;
    end
    FittedPhase = fminsearch(@RMSD_PhaseOffset,FittedPhase,[],Signal(FitStart:end));
    if nargin<2
        Phase=FittedPhase(1);
    else
        if isempty(Phase)
            Phase=FittedPhase(1);
        end
    end
    if sum(real(Signal*exp(1i*Phase)))<0
        Phase=Phase+pi;
    end
end

if FittedImaginaryOffset
    ImagOffset = FittedPhase(2);
else
    ImagOffset = 0;
end

%Do phase correction and normalize
if DirectDimension>1
    ReferenceSignal = (PrimaryData(1,:) - 1i*ImagOffset)*exp(1i*Phase);
    NormFactor = max(real(ReferenceSignal));
    sig = (PrimaryData(2,:) - 1i*ImagOffset)*exp(1i*Phase);
    correctedSignal = real(sig)./real(ReferenceSignal) + 1i*imag(ReferenceSignal)/NormFactor;
else
    correctedSignal = (PrimaryData - 1i*ImagOffset)*exp(1i*Phase);
    NormFactor = 1/max(real(correctedSignal));
    correctedSignal = NormFactor*correctedSignal;
end

end

function RMSD=RMSD_PhaseOffset(FitParam,ComplexSignal)
% Computes root mean square deviation of the imaginary part of
% phase-corrected data from zero before phase correction, an offset can be
% subtracted from the imaginary part
%
% FitParam(1)  phase correction phi (rad)
% FitParam(2)  offset
% ComplexSignal    complex data trace
%
% rmsd  root mean square deviation of imaginary part from zero
%
% G. Jeschke, 2009, Luis Fabregas 2019



if length(FitParam)>1
    if FitParam(2)<0
        RMSD = 1e6;
        return
    end
    %If requested, fit an imaginary offset additionally to phase
    ComplexSignal = ComplexSignal-1i*FitParam(2);
    ImaginaryComponent = imag(ComplexSignal*exp(1i*FitParam(1)));
    RMSD = norm(ImaginaryComponent);
    
else
    
    %Otherwise, just fit the correction phase
    ImaginaryComponent = imag(ComplexSignal*exp(1i*FitParam(1)));
    RMSD = norm(ImaginaryComponent);
    
end

end