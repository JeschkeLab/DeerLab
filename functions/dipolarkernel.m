% 
% DIPOLARKERNEL Computes the dipolar interaction kernel for the linear
%              transformation from distance-domain to time-domain
%
%       K = DIPOLARKERNEL(t,r) Computes the NxM point kernel for the 
%       trasnformation to the dipolar evolution function from the N-point 
%       time axis (t) and M-point distance axis (r).
%
%       K = DIPOLARKERNEL(t,r,B) Computes the kernel for the transformation
%       to the form factor function with a N-point background B multiplied.
%
%       K = DIPOLARKERNEL(T,R,[],'Property1',Value1,...) Computes the kernel 
%       for the transformation to the dipolar evolution function with the
%       input background B with additional optional arguments.
%
%       K = DIPOLARKERNEL(T,R,B,'Property1',Value1,...) Computes the kernel
%       for the transformation to the form factor function with a N-point
%       background B multiplied and additional optional arguments.
%
% Optional arguments can be specified by parameter/value pairs. The allowed
% properties to be passed as options can be set in any order. 
%
%   'KernelBType'   The way the background B is introduced into the kernel: 
%                      'full' - as passed without change
%                      'sqrt' - the sqrt of the background taken 
%                      'none' - background is not included in the kernel
%
%   'ExcitationBandwidth'   Excitation bandwith of the pulses in MHz to be
%                         used for limited bandwith excitation
%
%   'OvertoneCoeffs'    1D-Array of coefficients for correction of overtones
%                       in RIDME signals
%
%   'gValue'    g-value of the spin centers
%
%   'KernelCalcMethod'  The way the kernel is computed numerically:
%                       'fresnel' - uses Fresnel integrals for the kernel
%                       'explicit' - powder average computed explicitly
%
%   'Knots' Number knots for the grid of powder orientations to be used
%           in explicit kernel calculations    
%
%   'ModDepth' Modulation depth, if not specified it is determined from the
%              background fit
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.


function K = dipolarkernel(t,r,Background,ModDepth,varargin)

%--------------------------------------------------------------------------
%Input parsing
%--------------------------------------------------------------------------

if nargin<3 || isempty(Background)
    Background = ones(1,length(t));
end
if nargin<4 || isempty(ModDepth)
    ModDepth = 1;
end
if isa(Background,'char')
    varargin{end+1} = Background;
    varargin{end+1} = ModDepth;
    ModDepth = [];
    Background = ones(1,length(t));
end
%Check if user requested some options via name-value input
[KernelBType,ExcitationBandwidth,OvertoneCoeffs,gValue,KernelCalcMethod,Knots] = parseoptional({'KernelBType','ExcitationBandwidth','OvertoneCoeffs','gValue','KernelCalcMethod','Knots'},varargin);
%Validate the input variables
validKernelBTypes = ["none","full","sqrt"];
if isempty(KernelBType)
    KernelBType = 'full';
else
    validateattributes(KernelBType,{'char'},{},mfilename,'KernelBType')
    KernelBType = validatestring(KernelBType,validKernelBTypes);
end

if isempty(KernelCalcMethod)
    KernelCalcMethod = 'fresnel';
else
    validateattributes(KernelCalcMethod,{'char'},{'nonempty'},mfilename,'KernelCalcMethod')
end

if isempty(OvertoneCoeffs)
    OvertoneCoeffs = 1;
else
    validateattributes(OvertoneCoeffs,{'numeric'},{'nonnegative'},mfilename,'OvertoneCoeffs')
end

if ~isempty(ExcitationBandwidth)
    validateattributes(ExcitationBandwidth,{'numeric'},{'scalar','nonnegative'},mfilename,'ExcitationBandwidth')
end

if ~isempty(Background)
    validateattributes(Background,{'numeric'},{},mfilename,'Background')
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
if ~iscolumn(Background)
    Background = Background.';
end

if isempty(gValue)
    %Set g-value of isotropic nitroxide radical (old DA default)
    gValue = 2.004602204236924;
else
    validateattributes(gValue,{'numeric'},{'scalar','nonnegative'},mfilename,'gValue')
end

validateattributes(r,{'numeric'},{'nonempty','increasing','nonnegative'},mfilename,'DistanceAxis')
if numel(unique(round(diff(r),12)))~=1
    error('Distance axis must be a monotonically increasing vector.')
end
validateattributes(t,{'numeric'},{'nonempty','increasing'},mfilename,'TimeAxis')
checklengths(t,Background);

%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.LinkedHashMap;
end
hashKey = datahash({t,r,Background,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    K = java2mat(Output);
    return
end

%--------------------------------------------------------------------------
%Kernel construction
%--------------------------------------------------------------------------

Dimension = length(t);
TimeStep = mean(diff(t));

%Numerical dipolar frequency at 1 nm for given g-value, in MHz
nu0 = 51.92052556862238*gValue/2; % MHz nm^3
w0 = 2*pi*nu0; % Mrad s^-1 nm^3

%Convert time step to microseconds if given in nanoseconds
usesNanoseconds = TimeStep>=0.5;
if usesNanoseconds
    TimeStep = round(TimeStep)/1000;
    t = round(t)/1000;
end

%Get vector of dipolar frequencies
wdd = w0./(r.^3);

K = zeros(length(t),length(wdd));
for OvertoneIdx = 1:length(OvertoneCoeffs)
    
    switch KernelCalcMethod
        case 'explicit'
            % Method using explicit numerical powder average (slow)
            %----------------------------------------------------------
            %Pre-allocate cosine of powder angles
            costheta = linspace(0,1,Knots);
            %Sweep through all distances
            for DistanceIndex = 1:length(r)
                KernelTrace = 0;
                for theta = 1:Knots % average over cos(theta) angle (powder average)
                    KernelTrace = KernelTrace + cos(OvertoneIdx*wdd(DistanceIndex)*(1-3*costheta(theta)^2)*t);
                end
                % normalize dipolar time evolution trace
                K(:,DistanceIndex) = K(:,DistanceIndex) + OvertoneCoeffs(OvertoneIdx)*KernelTrace;
            end
            
        case  'fresnel'
            % Method using Fresnel integrals (fast)
            %----------------------------------------------------------
            %Compute dipolar kernel
            %Allocate products for speed
            wddt = OvertoneIdx*wdd.*t;
            kappa = sqrt(6*wddt/pi);
            %Compute Fresnel integrals of 0th order
            C = fresnelC(kappa);
            S = fresnelS(kappa);
            K = K + OvertoneCoeffs(OvertoneIdx)*sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S);
            %Replace undefined Fresnel NaN value at time zero
            K(isnan(K)) = 1;
    end
end

%If given, account for limited excitation bandwidth
if ~isempty(ExcitationBandwidth)
    K = exp(-wdd'.^2/ExcitationBandwidth^2).*K;
end

if strcmp(KernelCalcMethod,'explicit')
   [~,BckgStart] = min(abs(t));
    K = K./K(BckgStart,:);
end

% Build modulation depth and background into the kernel
%----------------------------------------------------------
if ModDepth~=1
  K = (1-ModDepth) + ModDepth*K;
end
switch KernelBType
    case 'none'
        % do not include B in K
    case 'full'
        K = K.*Background;
    case 'sqrt'
        K = K.*sqrt(Background);
end

%Normalize kernel
dr = mean(diff(r));
K = K*dr;

%Store output result in the cache
cachedData = addcache(cachedData,hashKey,K);

return