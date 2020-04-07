%
% DIPOLARSIGNAL Generate dipolar signal from distance distribution
%
%   V = DIPOLARSIGNAL(t,r,P)
%   Calculates the noiseless form factor (V) on a time axis (t) from the 
%   distance distribution (P) in a distance axis (r).
%
%   V = DIPOLARSIGNAL(t,r)
%   If no distribution is provided, then the dipolar signal corresponding 
%   to a single distance (r) is computed.
%
%   V = DIPOLARSIGNAL(...,'Property',Value)
%   Additional options can be passed as property-value pairs. You can specify
%   several name and value pair arguments in any order.
%
%   Name-Value Pair Arguments:
%
%   'ModDepth' - Modulation depth of the form factor
%
%   'B' - Array containing a background function
%
%   'NoiseLevel' - Level (standard deviation) of gaussian noise to add
%
%   'Scale' - Vertical scale to add to the ouput signal
%
%   'Phase' - Phase of the signal in radians.
%
%   'Overtones' - Array of RIDME overtone coefficients
%
%   'g' - Specifies the g-value of the electron spin center used to compute 
%              the dipolar frequencies from the given distance axis.
%
%   'Interference' - Cell array {A1 t1 A2 t2 ... @bg_bckg} containing the relative
%                    amplitudes and time shifts of the dipolar interferences. 
%                    The background model can be passed as a last argument to 
%                    include the time-shifted backgrounds
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [V] = dipolarsignal(t,r,P,varargin)

if nargin<3
    P = [];
    if length(r)~=1
       error('Only one distance is accepted if P is not provided.') 
    end
end
if ischar(P)
    varargin = [{P} varargin];
    P = [];
end
% Parse optional input arguments
[lambda,B,NoiseLevel,g,Scale,Overtones,Phase] = parseoptional({'ModDepth','Background','NoiseLevel','g','Scale','Overtones','Phase'},varargin);
% Validate inputs
if isempty(lambda)
    lambda = 1;
end
if isempty(B)
    B = ones(length(t),1);
end
if isempty(NoiseLevel)
    NoiseLevel = 0;
end
if isempty(Overtones)
    Overtones = [];
end
if isempty(Scale)
    Scale = 1;
end
if isempty(Phase)
    Phase = 0;
end
validateattributes(NoiseLevel,{'numeric'},{'scalar','nonnegative'},mfilename,'NoiseLevel')
validateattributes(lambda,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'ModDepth')
validateattributes(t,{'numeric'},{'nonempty'},mfilename,'t')
validateattributes(r,{'numeric'},{'increasing','nonempty','nonnegative'},mfilename,'r')
validateattributes(B,{'numeric'},{'2d'},mfilename,'B')
validateattributes(P,{'numeric'},{'2d'},mfilename,'P')
validateattributes(Scale,{'numeric'},{'scalar','nonnegative'},mfilename,'Offset')
validateattributes(Overtones,{'numeric'},{'2d','nonnegative'},mfilename,'Overtones')
validateattributes(Phase,{'numeric'},{'2d','scalar'},mfilename,'Phase')

if ~isempty(P) && length(r)~=length(P)
    error('The distance axis and distribution lengths must be equal.')
end

if lambda>1 || lambda<0
    error('Modulation depth must be in the range of 0 to 1.')
end
if numel(unique(round(diff(r),6)))~=1 && ~isscalar(r)
    error('Distance axis must be a monotonically increasing vector.')
end
if ~iscolumn(P)
   P = P.'; 
end

% Get the kernel
K = dipolarkernel(t,r,lambda,B,'OvertoneCoeffs',Overtones,'g',g);

% Calculate dipolar evolution function
if ~isempty(P)
    V = K*P;
else
    V = K;
end

% Phase rotate if given
if Phase~=0
    V = V.*exp(-1i*Phase);
end

% Add noise
if NoiseLevel>0
    V = V + whitegaussnoise(numel(t),NoiseLevel);
    if ~isreal(V)
        V = V + 1i*whitegaussnoise(numel(t),NoiseLevel);
    end
end

% Scale amplitude
V = V*Scale;

end