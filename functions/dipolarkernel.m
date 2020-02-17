%
% DIPOLARKERNEL Computes the dipolar interaction kernel for the linear
%              transformation from distance-domain to time-domain
%
%       K = DIPOLARKERNEL(t,r)
%       Computes the NxM point kernel for the trasnformation to the dipolar
%       evolution function from the N-point  time axis (t) in us or ns and M-point
%       distance axis (r) in ns or Angstrom.
%
%       K = DIPOLARKERNEL(t,r,lambda)
%       Computes the kernel for the transformation to the form factor function
%       with a modulation depth given by (lambda).
%
%       K = DIPOLARKERNEL(t,r,lambda,B)
%       Computes the kernel for the transformation to the form factor function
%       with a N-point background (B) multiplied.
%
%       K = DIPOLARKERNEL(t,r,'Property1',Value1,...)
%       K = DIPOLARKERNEL(t,r,lambda,'Property1',Value1,...)
%       K = DIPOLARKERNEL(t,r,lambda,B,'Property1',Value1,...)
%       Optional arguments can be specified by parameter/value pairs.
%
% The following properties are available.
%
%   'ExcitationBandwidth' - Excitation bandwith of the pulses in MHz to be
%                           used for limited bandwith excitation
%
%   'OvertoneCoeffs' - 1D-Array of coefficients for correction of overtones
%                      in RIDME signals
%
%   'g' - g-value of the spin centers
%
%   'Method' - The way the kernel is computed numerically:
%               'fresnel' - uses Fresnel integrals for the kernel (default)
%               'grid' - powder average computed explicitly (slow)
%
%   'nKnots' - Number knots for the grid of powder orientations to be used
 %             in explicit kernel calculations
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function K = dipolarkernel(t,r,varargin)

% Input parsing
%-------------------------------------------------------------------------------
if nargin<2
    error('At least two inputs (t and r) are required.');
end

% Set default parameters
lambdaT0 = 1;
B = [];
proplist = varargin;

if numel(proplist)>=1 && ~ischar(proplist{1})
    lambdaT0 = proplist{1};
    proplist(1) = [];
end

if numel(proplist)>=1 && ~ischar(proplist{1})
    B = proplist{1};
    proplist(1) = [];
end

% Check if user requested some options via name-value input
[ExcitationBandwidth,OvertoneCoeffs,g,Method,nKnots,useCache] = ...
    parseoptional({'ExcitationBandwidth','OvertoneCoeffs','g','Method','nKnots','Cache'},proplist);
if isempty(useCache)
    useCache = true;
end
if isempty(Method)
    Method = 'fresnel';
end
validateattributes(Method,{'char'},{'nonempty'},mfilename,'Method');

if isempty(OvertoneCoeffs)
    OvertoneCoeffs = 1;
end
validateattributes(OvertoneCoeffs,{'numeric'},{'nonnegative'},mfilename,'OvertoneCoeffs');

if ~isempty(ExcitationBandwidth)
    validateattributes(ExcitationBandwidth,{'numeric'},{'scalar','nonnegative'},mfilename,'ExcitationBandwidth');
end

if ~isempty(B)
    validateattributes(B,{'numeric','function_handle'},{},mfilename,'B');
    if isnumeric(B) 
        checklengths(t,B);
    end
end

if isempty(nKnots)
    nKnots = 5001;
end
validateattributes(nKnots,{'numeric'},{'scalar'},mfilename,'Knots');

% Make sure all vectors are column vectors
t = t(:);
r = r(:);
if isnumeric(B)
    B = B(:);
end

ge = 2.00231930436256; % free-electron g factor (CODATA 2018 value)
if isempty(g)
    g = ge;
end
validateattributes(g,{'numeric'},{'scalar','nonnegative'},mfilename,'g');

validateattributes(r,{'numeric'},{'nonempty','increasing','nonnegative'},mfilename,'r');
if numel(unique(round(diff(r),6)))~=1 && length(r)~=1
    error('Distance axis must be a monotonically increasing vector.');
end
validateattributes(t,{'numeric'},{'nonempty'},mfilename,'t');

% Memoization
%-------------------------------------------------------------------------------
persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({t,r,varargin});
if cachedData.containsKey(hashKey)  && useCache
    Output = cachedData.get(hashKey);
    K = java2mat(Output);
    return
end

% Validation of the multipathway parameters
%-------------------------------------------------------------------------------
if numel(lambdaT0)==1
  lambda = lambdaT0;
  T0 = 0;
else
  lambda = lambdaT0(:,1);
  T0 = lambdaT0(:,2);
end

% Assert that amplitude of non-modulated pathways is not negative
if sum(lambda)>1
    error('Sum of all lambdas cannot be larger than 1.');
end
lambda0 = 1 - sum(lambda);

nPathways = numel(lambda);
if nPathways>1
    if ~isempty(B) && ~isa(B,'function_handle')
        error('For a multi-pathway model, B must be a function handle.');
    end
end

% Kernel construction
%-------------------------------------------------------------------------------
% Numerical dipolar frequency at 1 nm for given g-value, in MHz
muB = 9.2740100783e-24; % Bohr magneton, J/T (CODATA 2018 value);
mu0 = 1.25663706212e-6; % magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
h = 6.62607015e-34; % Planck constant, J/Hz (CODATA 2018)
nu0 = (mu0/4/pi)*(muB*g)^2/h*1e21; % MHz nm^3
w0 = 2*pi*nu0; % Mrad s^-1 nm^3

% Get vector of dipolar frequencies at all distances
wdd = w0./r.^3;

% Set kernel matrix calculation method
switch Method
    case 'fresnel'
        kernelmatrix = @(t)kernelmatrix_fresnel(t,wdd,OvertoneCoeffs);
    case 'grid'
        kernelmatrix = @(t)kernelmatrix_grid(t,wdd,OvertoneCoeffs,nKnots);
end

% Build dipolar signal, summing over all pathways
nPathways = numel(lambda);
K = lambda0;
for p = 1:nPathways
    K = K + lambda(p)*kernelmatrix(t-T0(p));
end

% Multiply by background
if isa(B,'function_handle')
    for p = 1:nPathways
        K = K.*B(t-T0(p),lambda(p));
    end
else
    if ~isempty(B)
        K = K.*B;
    end
end

% Include dr into kernel
if length(r)>1
    dr = mean(diff(r));
    K = K*dr;
end

% If given, include limited excitation bandwidth
if ~isempty(ExcitationBandwidth)
    K = exp(-wdd.'.^2/ExcitationBandwidth^2).*K;
end

% Store output result in the cache
cachedData = addcache(cachedData,hashKey,K);

end

%===============================================================================
% Calculate kernel using Fresnel integrals (fast)
function K = kernelmatrix_fresnel(t,wdd,OvertoneCoeffs)

K = zeros(length(t),length(wdd));

for n = 1:numel(OvertoneCoeffs)
    ph = n*wdd.'.*abs(t);
    kappa = sqrt(6*ph/pi);
    C = fresnelC(kappa)./kappa;
    S = fresnelS(kappa)./kappa;
    K = K + OvertoneCoeffs(n)*(cos(ph).*C + sin(ph).*S);
end

% Replace NaN values at time zero with 1
K(isnan(K)) = 1;

end

%===============================================================================
% Calculate kernel using grid-based powder integration (slow)
% (converges very slowly with nKnots)
function K = kernelmatrix_grid(t,wdd,OvertoneCoeffs,nKnots)

K = zeros(length(t),length(wdd));

costheta = linspace(0,1,nKnots);
q = 1 - 3*costheta.^2;

for n = 1:numel(OvertoneCoeffs)
    
    for ir = 1:numel(wdd)
        K_ = 0;
        for itheta = 1:nKnots
            K_ = K_ + cos(n*wdd(ir)*q(itheta)*t);
        end
        K(:,ir) = K(:,ir) + OvertoneCoeffs(n)*K_/nKnots;
    end
    
end

end
