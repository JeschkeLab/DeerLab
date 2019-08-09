function [Kernel] = dipolarkernel(TimeAxis,DistanceAxis,Background,varargin)

%--------------------------------------------------------------------------
%Input parsing
%--------------------------------------------------------------------------
%Check if user requested some options via name-value input
[KernelBType,ExcitationBandwidth,OvertoneCoeffs,gValue,KernelCalcMethod,Knots] = parseoptional({'KernelBType','ExcitationBandwidth','OvertoneCoeffs','gValue','KernelCalcMethod','Knots'},varargin);
%Validate the input variables
if nargin<3 || isempty(Background)
    Background = ones(1,length(TimeAxis));
end

validKernelBTypes = ["none","full","sqrt"];
if isempty(KernelBType)
    KernelBType = 'sqrt';
else
    validateattributes(KernelBType,{'char'},{},mfilename,'KernelBType')
    validatestring(KernelBType,validKernelBTypes);
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

if ~iscolumn(TimeAxis)
    TimeAxis = TimeAxis.';
end

if iscolumn(DistanceAxis)
    DistanceAxis = DistanceAxis.';
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

validateattributes(DistanceAxis,{'numeric'},{'nonempty','increasing','nonnegative'},mfilename,'DistanceAxis')
validateattributes(TimeAxis,{'numeric'},{'nonempty','increasing','nonnegative'},mfilename,'TimeAxis')
checklengths(TimeAxis,Background);

%--------------------------------------------------------------------------
%Memoization
%--------------------------------------------------------------------------

persistent cachedData
if isempty(cachedData)
    cachedData =  java.util.Hashtable;
end
hashKey = datahash({TimeAxis,DistanceAxis,Background,varargin});
if cachedData.containsKey(hashKey)
    Output = cachedData.get(hashKey);
    Kernel = java2mat(Output);
    return
end

%--------------------------------------------------------------------------
%Kernel construction
%--------------------------------------------------------------------------

Dimension = length(TimeAxis);
TimeStep = mean(diff(TimeAxis));

%Numerical dipolar frequency at 1 nm for given g-value
ny0 = 51.92052556862238*gValue/2;

%Convert time step to miliseconds if given in nanosecond
if TimeStep>1
    TimeStep = round(TimeStep)/1000;
    TimeAxis = round(TimeAxis)/1000;
end

%Numerical angular dipolar frequency at 1 nm for g=ge
w0 = 2*pi*ny0;
%Get vector of dipolar frequencies
wdd = w0./(DistanceAxis.^3);

Kernel = zeros(length(TimeAxis),length(wdd));
for OvertoneIdx=1:length(OvertoneCoeffs)
    
    switch KernelCalcMethod
        case 'explicit'
            % Method using explicit numerical powder average (slow)
            %----------------------------------------------------------
            %Pre-allocate cosine of powder angles
            costheta = linspace(0,1,Knots);
            %Sweep through all distances 
            for DistanceIndex = 1:length(DistanceAxis)
                KernelTrace = 0;
                for theta = 1:Knots % average over cos(theta) angle (powder average)
                    KernelTrace = KernelTrace + cos(OvertoneIdx*wdd(DistanceIndex)*(1-3*costheta(theta)^2)*TimeAxis);
                end
                % normalize dipolar time evolution trace
                Kernel(:,DistanceIndex) = Kernel(:,DistanceIndex) + OvertoneCoeffs(OvertoneIdx)*KernelTrace;
            end
        case  'fresnel'
            % Method using Fresnel integrals (fast)
            %----------------------------------------------------------
            %Compute dipolar kernel
            %Allocate products for speed
            wddt = OvertoneIdx*wdd.*TimeAxis;
            kappa = sqrt(6*wddt/pi);
            %Compute Fresnel integrals of 0th order
            C = fresnelC(kappa);
            S = fresnelS(kappa);
            Kernel = Kernel + OvertoneCoeffs(OvertoneIdx)*sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S);
            %Replace undefined Fresnel NaN value at time zero
            Kernel(1,:) = 1;
    end
end

%Normalize kernel
Kernel = Kernel./Kernel(1,:);

%If given, account for limited excitation bandwidth
if ~isempty(ExcitationBandwidth)
    Kernel = exp(-wdd'.^2/ExcitationBandwidth^2).*Kernel;
end

% Build the background into the kernel
%----------------------------------------------------------
if ~all(Background == 1)
    ModDepth = 1/Background(1) - 1;
    Kernel = ModDepth*Kernel + 1;
end
switch KernelBType
    case 'none'
        Background = ones(Dimension,1);
    case 'sqrt'
        Background = sqrt(Background);
end
Kernel = Kernel.*Background;

%Store output result in the cache
cachedData = addcache(cachedData,hashKey,Kernel);

return