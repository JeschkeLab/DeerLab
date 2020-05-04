%
% DIPOLARSIGNAL Generate dipolar signal from distance distribution
%
%    V = DIPOLARSIGNAL(t,r)
%    V = DIPOLARSIGNAL(t,r,P)
%    V = DIPOLARSIGNAL(t,r,P,lambda)
%    V = DIPOLARSIGNAL(t,r,P,lambda,B)
%    V = DIPOLARSIGNAL(t,r,P,pathinfo)
%    V = DIPOLARSIGNAL(t,r,P,pathinfo,B)
%    V = DIPOLARSIGNAL(___,'Property',value,___)
%
%   Calculates the dipolar time-domain signal (V) on a time axis (t), in
%   microseconds, from the distance distribution (P) on a distance axis (r),
%   in nanometers. If given, the modulation amplitude (lambda) and a
%   background (B) are included. Multiple pathways can be specified in
%   (pathinfo). Additional options can be given by name-value pairs.
%
%   If no distribution is provided, then the dipolar signal corresponding 
%   to a single distance (r) is computed.
%
%   Name-value pair arguments:
%
%   'ModDepth'   - Modulation depth of the form factor
%   'Background' - Array containing a background function
%   'NoiseLevel' - Level (standard deviation) of gaussian noise to add
%   'Scale'      - Vertical scale to add to the ouput signal
%   'Phase'      - Phase of the signal in radians.
%   'Overtones'  - Array of RIDME overtone coefficients
%   'g'          - Specifies the g-value of the electron spin center used to
%                  compute the dipolar frequencies from the given distance axis.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function V = dipolarsignal(t,r,varargin)

if nargin<2
  error('At least two inputs (t and r) are required.');
end

% Set default parameters
proplist = varargin;
P = [];
pathinfo = [];
B = [];

if numel(proplist)>=1 && ~ischar(proplist{1})
    P = proplist{1};
    proplist(1) = [];
end
if numel(proplist)>=1 && ~ischar(proplist{1})
    pathinfo = proplist{1};
    proplist(1) = [];
end
if numel(proplist)>=1 && ~ischar(proplist{1})
    B = proplist{1};
    proplist(1) = [];
end

% Parse optional input arguments
names = {'NoiseLevel','g','Scale','Overtones','Phase'};
[NoiseLevel,g,Scale,Overtones,Phase] = parseoptional(names,proplist);

% Validate inputs
if isempty(pathinfo)
    pathinfo = 1;
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
%validateattributes(lambda,{'numeric'},{'scalar','nonnegative','nonempty'},mfilename,'pathinfo')
if ~isa(B,'function_handle')
    validateattributes(B,{'numeric'},{'2d'},mfilename,'B')
end
validateattributes(NoiseLevel,{'numeric'},{'scalar','nonnegative'},mfilename,'NoiseLevel')
validateattributes(t,{'numeric'},{'nonempty'},mfilename,'t')
validateattributes(r,{'numeric'},{'increasing','nonempty','nonnegative'},mfilename,'r')
validateattributes(P,{'numeric'},{'2d'},mfilename,'P')
validateattributes(Scale,{'numeric'},{'scalar','nonnegative'},mfilename,'Offset')
validateattributes(Overtones,{'numeric'},{'2d','nonnegative'},mfilename,'Overtones')
validateattributes(Phase,{'numeric'},{'2d','scalar'},mfilename,'Phase')

if ~isempty(P) && numel(r)~=numel(P)
    error('The distance axis and distribution lengths must be equal.')
end

if numel(unique(round(diff(r),6)))~=1 && ~isscalar(r)
    error('Distance axis must be a monotonically increasing vector.')
end
if ~iscolumn(P)
   P = P.';
end

% Assemble pathway info
if numel(pathinfo)==1
    lambda = pathinfo;
    if numel(lambda)~=1 || lambda>1 || lambda<0
      error('Modulation depth must be a number in the range of 0 to 1.')
    end
    pathinfo = [lambda 0; 1-lambda NaN];
end

% Get the kernel matrix
K = dipolarkernel(t,r,pathinfo,B,'OvertoneCoeffs',Overtones,'g',g);

% Calculate dipolar evolution function
if ~isempty(P)
    V = K*P;
else
    V = K;
end

% Rotate phase, add noise, scale
%-------------------------------------------------------------------------------
% Phase rotate
if Phase~=0
    V = V.*exp(-1i*Phase);
end

% Add noise
if NoiseLevel~=0
    V = V + whitegaussnoise(numel(t),NoiseLevel);
    if ~isreal(V)
        V = V + 1i*whitegaussnoise(numel(t),NoiseLevel);
    end
end

% Scale amplitude
if Scale~=1
    V = V*Scale;
end

end
