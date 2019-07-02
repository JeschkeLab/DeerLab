function [correctedSignal,Phase,ImagOffset] = correctPhase(PrimaryData,Phase,FittedImaginaryOffset)

if nargin<3
  FittedImaginaryOffset = false;
end

if nargin==1 || nargin==3
    Phase = [];
end
DirectDimension = size(PrimaryData,1);

if DirectDimension>1
  Signal=PrimaryData(1,:);
else
  Signal=PrimaryData;
end

% If phse is not provided, then fit it
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

% make phase correction and normalize
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