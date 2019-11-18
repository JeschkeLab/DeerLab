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
% The allowed properties to be passed as options can be set in any order.
%
%   'ExcitationBandwidth' - Excitation bandwith of the pulses in MHz to be
%                           used for limited bandwith excitation
%
%   'OvertoneCoeffs' - 1D-Array of coefficients for correction of overtones
%                      in RIDME signals
%
%   'gValue' - g-value of the spin centers
%
%   'Method' - The way the kernel is computed numerically:
%                       'fresnel' - uses Fresnel integrals for the kernel
%                       'explicit' - powder average computed explicitly
%
%   'Knots' - Number knots for the grid of powder orientations to be used
%             in explicit kernel calculations
%
%   'Interference' - Cell array {A1 t1 A2 t2 ... @td_bckg} containing the relative
%                    amplitudes and time shifts of the dipolar interferences. 
%                    The background model can be passed as a last argument to 
%                    include the time-shifted backgrounds
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.



function K = dipolarkernel(t,r,lambda,B,varargin)

% Input parsing
%--------------------------------------------------------------------------
switch nargin
    case {1,2}
        lambda = [];
        B = [];
end

% Case ModDepth and B not passed
if nargin>2 && isa(lambda,'char')
    varargin = [{lambda} {B} varargin];
    lambda = [];
    B = [];
end
if nargin < 4
    B = [];
end
% Case ModDepth passed but not B
if nargin>3 && isa(B,'char')
    varargin = [{B} varargin];
    B = [];
end

% Check if user requested some options via name-value input
[ExcitationBandwidth,OvertoneCoeffs,gValue,Method,Knots,InterferenceParam] = ...
    parseoptional({'ExcitationBandwidth','OvertoneCoeffs','gValue','Method','Knots','Interference'},varargin);

if isempty(Method)
    Method = 'fresnel';
else
    validateattributes(Method,{'char'},{'nonempty'},mfilename,'Method');
end

if isempty(OvertoneCoeffs)
    OvertoneCoeffs = 1;
else
    validateattributes(OvertoneCoeffs,{'numeric'},{'nonnegative'},mfilename,'OvertoneCoeffs');
end

if ~isempty(ExcitationBandwidth)
    validateattributes(ExcitationBandwidth,{'numeric'},{'scalar','nonnegative'},mfilename,'ExcitationBandwidth');
end

if ~isempty(B)
    validateattributes(B,{'numeric'},{},mfilename,'B')
    checklengths(t,B);
end

if isempty(Knots)
    Knots = 1001;
else
    validateattributes(Knots,{'numeric'},{'scalar'},mfilename,'Knots')
end

if ~iscolumn(t)
    t = t.';
end

if iscolumn(r)
    r = r.';
end
if ~iscolumn(B)
    B = B.';
end

if isempty(gValue)
    %Set g-value of isotropic nitroxide radical (old DA default)
    gValue = 2.004602204236924;
else
    validateattributes(gValue,{'numeric'},{'scalar','nonnegative'},mfilename,'gValue')
end

validateattributes(r,{'numeric'},{'nonempty','increasing','nonnegative'},mfilename,'r')
if numel(unique(round(diff(r),6)))~=1 && length(r)~=1
    error('Distance axis must be a monotonically increasing vector.')
end
validateattributes(t,{'numeric'},{'nonempty','increasing'},mfilename,'t')

% Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({t,r,lambda,B,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    K = java2mat(Output);
    return
end

% Kernel construction
%--------------------------------------------------------------------------

dt = mean(diff(t));

% Numerical dipolar frequency at 1 nm for given g-value, in MHz
nu0 = 51.92052556862238*gValue/2; % MHz nm^3
w0 = 2*pi*nu0; % Mrad s^-1 nm^3

% Convert time step to microseconds if given in nanoseconds
if isnanosecond(t)
    t = round(t,12)/1000; % ns->us
end

% Convert distance axis to nanoseconds if givne in Angstrom
rinput = r;
if ~isnanometer(r)
    r = r/10;
end
if length(r)~=1
    dr = mean(diff(rinput));
else
    dr = 1;
end

if ~isempty(InterferenceParam)
    if ~iscell(InterferenceParam)
        InterferenceParam = num2cell(InterferenceParam);
    end
    if numel(InterferenceParam)<2
        error('Interference option must contain at least two elements [amp tau]')
    end
    if mod(numel(InterferenceParam),2)~=0 && ~isa(InterferenceParam{end},'function_handle')
        error('The last element of the Interference option must be valid function handle.')
    elseif mod(numel(InterferenceParam),2)~=0 && isa(InterferenceParam{end},'function_handle')
        Bmodel = InterferenceParam{end};
        InterferenceParam(end) = [];
    else
        Bmodel = [];
    end
end
% Get absolute time axis scale (required for negative times)
traw = t;
t = abs(t); %ns

% Get vector of dipolar frequencies
wdd = w0./(r.^3);

K = zeros(length(t),length(wdd));
for OvertoneIdx = 1:numel(OvertoneCoeffs)
    
    switch Method
        case 'explicit'
            % Method using explicit numerical powder average (slow)
            %----------------------------------------------------------
            % Pre-allocate cosine of powder angles
            costheta = linspace(0,1,Knots);
            % Sweep through all distances
            for DistanceIndex = 1:length(r)
                KTrace = 0;
                for theta = 1:Knots % average over cos(theta) angle (powder average)
                    KTrace = KTrace + cos(OvertoneIdx*wdd(DistanceIndex)*(1-3*costheta(theta)^2)*t);
                end
                % normalize dipolar time evolution trace
                K(:,DistanceIndex) = K(:,DistanceIndex) + OvertoneCoeffs(OvertoneIdx)*KTrace;
            end
            
        case  'fresnel'
            % Method using Fresnel integrals (fast)
            %----------------------------------------------------------
            % Compute dipolar kernel
            % Allocate products for speed
            wddt = OvertoneIdx*wdd.*t;
            kappa = sqrt(6*wddt/pi);
            % Compute Fresnel integrals of 0th order
            C = fresnelC(kappa);
            S = fresnelS(kappa);
            K = K + OvertoneCoeffs(OvertoneIdx)*sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S);
            % Replace undefined Fresnel NaN value at time zero
            K(isnan(K)) = 1;
    end
end

% If given, account for limited excitation bandwidth
if ~isempty(ExcitationBandwidth)
    K = exp(-wdd'.^2/ExcitationBandwidth^2).*K;
end

if strcmp(Method,'explicit')
    [~,BckgStart] = min(t);
    K = K./K(BckgStart,:);
end

% Build modulation depth and background into the kernel
%----------------------------------------------------------
if ~isempty(lambda)
    K = (1-lambda) + lambda*K;
end
if ~isempty(B)
    K = K.*B;
end

% Build dipolar interference model
%----------------------------------------------------------
Kinter = zeros(length(t),length(wdd));
if ~isempty(InterferenceParam)
    %Loop over all interferences
    for i=1:2:numel(InterferenceParam)
        %Extract current interference parameters
        amp = InterferenceParam{i};
        tau = InterferenceParam{i + 1};
        %If background model is given then compute the background
        if ~isempty(Bmodel)
            %Get background parameters (non-elegant but efficient)
            [~,~,Bparam] = fitbackground(B,traw,Bmodel,[min(t) max(t)]);
            %Compute time-shifted background
            Bshifted = Bmodel(abs(traw - tau),Bparam);
        else
            Bshifted = ones(size(t));
        end
        %Add the scaled dipolar kernel for the current shifted time axis
        Kinter = Kinter + amp*dipolarkernel(traw - tau,r,lambda,Bshifted)/dr;
    end
end
%Add interference kernels and correct for understimated modulation depth
K = K + Kinter - Kinter(1,:);

% Normalize kernel
K = K*dr;

%Store output result in the cache
cachedData = addcache(cachedData,hashKey,K);

return
