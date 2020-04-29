%
% DIPOLARKERNEL Dipolar kernel matrix
%
%   K = DIPOLARKERNEL(t,r)
%   K = DIPOLARKERNEL(t,r,lambda)
%   K = DIPOLARKERNEL(t,r,lambda,B)
%   K = DIPOLARKERNEL(t,r,pathinfo)
%   K = DIPOLARKERNEL(t,r,pathinfo,B)
%   K = DIPOLARKERNEL(___,'Property',value,___)
%
%  Computes the NxM kernel matrix K for the transformation of a distance
%  distribution to a dipolar evolution function, using the M-point distance
%  axis r (in nanometers) and the N-point time axis t (in microseconds).
%  If the modulation depth lambda, is given, it is included in K.
%  If an N-point background (B) is given, it is included in K as well.
%  Optional arguments can be specified by name-value pairs.
%
%  Inputs:
%     t         N-element time axis, in microseconds
%     r         M-element distance axis, in nanometers
%     lambda    modulation depth (between 0 and 1)
%               this is equivalent to pathinfo = [1-lambda NaN 0; lambda 0 1]
%     pathinfo  px2 or px3 array of modulation depths lambda, refocusing points
%               T0, and harmonics n for multiple pathways, each row contains
%               [lambda T0 n] or [lambda T0] for one pathway. If n is not given
%               it is assumed to be 1.
%     B         N-element array with background decay, or a function handle
%               for a background model: @(t,lambda)bg_model(t,par,lambda)
%
%  Outputs:
%     K         NxM kernel matrix, such that the time-domain signal for a
%               Mx1 distance distribution vector P is K*P
%
%  Name-value pairs:
%
%   'ExcitationBandwidth' - Excitation bandwith of the pulses in MHz to be
%                           used for limited bandwith excitation
%   'OvertoneCoeffs' - 1D array of coefficients for overtones in RIDME signals
%   'g'              - g-values of the spin centers, [g1 g2]
%   'Method'         - Numerical method for kernel matrix calculation:
%                      'fresnel' - uses Fresnel integrals for the kernel (default)
%                      'integral' - uses MATLAB's integral() function (slow, accurate)
%                      'grid' - powder average computed explicitly (slow, inaccurate)
%   'nKnots'         - Number of knots for the grid of powder orientations to be used
%                      in the 'grid' kernel calculation method
%   'Renormalize'    - Re-normalization of multi-pathway kernels to ensure the
%                      equality V(0) == 1 is satisfied. (default = true)
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function K = dipolarkernel(t,r,varargin)

% Input parsing
%-------------------------------------------------------------------------------
if nargin<2
    error('At least two inputs (t and r) are required.');
end

% Set default parameters
pathinfo = 1;
B = [];
proplist = varargin;

if numel(proplist)>=1 && ~ischar(proplist{1})
    if ~isempty( proplist{1})
        pathinfo = proplist{1};
    end
    proplist(1) = [];
end

if numel(proplist)>=1 && ~ischar(proplist{1})
    B = proplist{1};
    proplist(1) = [];
end

% Check if user requested some options via name-value input
[ExcitationBandwidth,OvertoneCoeffs,g,Method,nKnots,useCache,Renormalize] = ...
    parseoptional({'ExcitationBandwidth','OvertoneCoeffs','g','Method','nKnots','Cache','Renormalize'},proplist);
if isempty(useCache)
    useCache = true;
end
if isempty(Method)
    Method = 'fresnel';
end
validateattributes(Method,{'char'},{'nonempty'},mfilename,'Method');

if isempty(Renormalize)
    Renormalize = true;
else
    validateattributes(Renormalize,{'logical'},{'nonempty'},mfilename,'Renormalize');
end
if isempty(OvertoneCoeffs)
    OvertoneCoeffs = 1;
end
validateattributes(OvertoneCoeffs,{'numeric'},{'nonnegative'},mfilename,'OvertoneCoeffs');

if ~isempty(ExcitationBandwidth)
    validateattributes(ExcitationBandwidth,{'numeric'},{'scalar','nonnegative'},mfilename,'ExcitationBandwidth');
    if strcmp(Method,'fresnel')
       error('Compensation for limited excitation bandwidth is not compatible with the ''fresnel'' method. Please use the ''integral'' or '' grid'' methods.') 
    end
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
    g = [ge ge];
end
validateattributes(g,{'numeric'},{'nonempty','nonnegative'},mfilename,'g');
if numel(g)~=2
    error('The array supplied for ''g'' must contain one or two elements.');
end

validateattributes(r,{'numeric'},{'nonempty','increasing','nonnegative'},mfilename,'r');
if numel(unique(round(diff(r),6)))~=1 && length(r)~=1
    error('Distance axis must be a monotonically increasing vector.');
end
validateattributes(t,{'numeric'},{'nonempty'},mfilename,'t');


% Validation of the multi-pathway parameters
%-------------------------------------------------------------------------------
if ~isnumeric(pathinfo) || ~isreal(pathinfo)
    error('lambda/pathinfo must be a numeric array.');
end

if numel(pathinfo)==1
    lambda = pathinfo;
    pathinfo = [1-lambda NaN; lambda 0];
end

if ~any(size(pathinfo,2)==[2 3])
  error('pathinfo must be a numeric array with two or three columns.');
end
if any(isnan(pathinfo(:,1)))
  error('In pathinfo, NaN can only appear in the second column (refocusing time) e.g. path(1,:) = [Lam0 NaN];');
end

% Normalize the pathway amplitudes to unity
pathinfo(:,1) = pathinfo(:,1)/sum(pathinfo(:,1));
lambda = pathinfo(:,1);
T0 = pathinfo(:,2);
if size(pathinfo,2)==2
  n = ones(size(T0));
else
  n = pathinfo(:,3);
end

% Combine all unmodulated components into Lambda0, and eliminate from list
unmodulated = isnan(T0);
Lambda0 = sum(lambda(unmodulated));
lambda(unmodulated) = [];
T0(unmodulated)  = [];
n(unmodulated) = [];

% Fold overtones into pathway list
nCoeffs = numel(OvertoneCoeffs);
if nCoeffs>0
  lambda = reshape(lambda*OvertoneCoeffs(:).',[],1);
  T0 = reshape(repmat(T0,1,nCoeffs),[],1);
  n = reshape(n*(1:nCoeffs),[],1);
end

nModPathways = numel(lambda);
if nModPathways>1
    if ~isempty(B) && (~isa(B,'function_handle') || nargin(B)~=2)
        error(['For a model with multiple modulated pathways, B must be a ',...
               'function handle of the type: @(t,lambda) bg_model(t,par,lambda)']);
    end
end

% Kernel matrix construction
%-------------------------------------------------------------------------------
% Numerical dipolar frequency at 1 nm for given g-value, in MHz
muB = 9.2740100783e-24; % Bohr magneton, J/T (CODATA 2018 value);
mu0 = 1.25663706212e-6; % magnetic constant, N A^-2 = T^2 m^3 J^-1 (CODATA 2018)
h = 6.62607015e-34; % Planck constant, J/Hz (CODATA 2018)
nu0 = (mu0/4/pi)*muB^2*g(1)*g(2)/h*1e21; % MHz nm^3
w0 = 2*pi*nu0; % Mrad s^-1 nm^3

% Get vector of dipolar frequencies at all distances
wdd = w0./r.^3;

% Memoization
if useCache
    kernelmatrix_ = memoize(@calckernelmatrix);
else
    kernelmatrix_ = @calckernelmatrix;
end
kernelmatrix = @(t)kernelmatrix_(Method,t,wdd,ExcitationBandwidth,nKnots);

% Build dipolar kernel matrix, summing over all pathways
K = Lambda0;
Knorm = Lambda0;
for p = 1:nModPathways
    K = K + lambda(p)*kernelmatrix(n(p)*(t-T0(p)));
    Knorm = Knorm + lambda(p)*kernelmatrix(-T0(p)*n(p));
end
if Renormalize
    K = K./Knorm;
end

% Multiply by background(s)
if isa(B,'function_handle')
    B_ = dipolarbackground(t,pathinfo,B,'OvertoneCoeffs',OvertoneCoeffs,'Renormalize',Renormalize);
    K = K.*B_;
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

end


%===============================================================================
% Calculate elementary kernel
function K = calckernelmatrix(method,t,wdd,wex,nKnots)

switch method
    case 'fresnel'
        K = kernelmatrix_fresnel(t,wdd);
    case 'grid'
        K = kernelmatrix_grid(t,wdd,wex,nKnots);
    case 'integral'
        K = kernelmatrix_integral(t,wdd,wex);
    otherwise
        error('Kernel calculation method ''%s'' is unknown.',method);
end

end

%===============================================================================
% Calculate kernel using Fresnel integrals (fast)
function K = kernelmatrix_fresnel(t,wdd)

ph = wdd.'.*abs(t);
kappa = sqrt(6*ph/pi);
C = fresnelC(kappa)./kappa;
S = fresnelS(kappa)./kappa;
K = cos(ph).*C + sin(ph).*S;

% Replace NaN values at time zero with 1
K(isnan(K)) = 1;

end

%===============================================================================
% Calculate kernel using grid-based powder integration (slow)
% (converges very slowly with nKnots)
function K = kernelmatrix_grid(t,wdd,wex,nKnots)

K = zeros(numel(t),numel(wdd));

costheta = linspace(0,1,nKnots);

q = 1-3*costheta.^2;
for ir = 1:numel(wdd)
  D_ = 0;
  for itheta = 1:nKnots
    C = cos(wdd(ir)*q(itheta)*abs(t));
    % If given, include limited excitation bandwidth
    if ~isempty(wex)
        C = C*exp(-(wdd(ir)*q(itheta)).^2/wex^2);
    end
     D_ = D_ + C;
  end
  K(:,ir) = D_/nKnots;
end

end

%===============================================================================
% Calculate kernel using MATLAB's integrator (accurate but slow)
function K = kernelmatrix_integral(t,wdd,wex)

K = zeros(numel(t),numel(wdd));

for ir = 1:numel(wdd)
    % If given, include limited excitation bandwidth
    if ~isempty(wex)
        fun = @(costheta) cos(wdd(ir)*abs(t)*(1-3*costheta.^2))*exp(-(wdd(ir)*(1-3*costheta.^2))^2/wex^2);
    else
        fun = @(costheta) cos(wdd(ir)*abs(t)*(1-3*costheta.^2));
    end
    K(:,ir) = integral(fun,0,1,'ArrayValued',true,'AbsTol',1e-6,'RelTol',1e-6);
end

end
