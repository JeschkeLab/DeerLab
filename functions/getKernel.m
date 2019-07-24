function [Kernel] = getKernel(TimeAxis,DistanceAxis,Background,varargin)

%--------------------------------------------------------------------------
%Input parsing
%--------------------------------------------------------------------------
%Check if user requested some options via name-value input
[KernelBType,ExcitationBandwidth,OvertoneCoeffs,gValue] = parseOptional({'KernelBType','ExcitationBandwidth','OvertoneCoeffs','gValue'},varargin);
%Validate the input variables
if nargin<3 || isempty(Background)
    Background = ones(1,length(TimeAxis));
end
validKernelBTypes = ["none","full","sqrt"];
if isempty(KernelBType)
    KernelBType = 'sqrt';
end
if isempty(OvertoneCoeffs)
    OvertoneCoeffs = 1;
else
    validateattributes(OvertoneCoeffs,{'numeric'},{'nonnegative'},mfilename,'OvertoneCoeffs')
end
if ~isempty(ExcitationBandwidth)
    validateattributes(ExcitationBandwidth,{'numeric'},{'scalar','nonnegative'},'getKernel','ExcitationBandwidth')
end
if ~isempty(Background)
    validateattributes(Background,{'numeric'},{},'getKernel','Background')
end
if iscolumn(TimeAxis)
    TimeAxis = TimeAxis';
end
if ~iscolumn(DistanceAxis)
    DistanceAxis = DistanceAxis';
end
if iscolumn(Background)
    Background = Background';
end
if isempty(gValue)
    %Set g-value of isotropic nitroxide radical (old DA default)
    gValue = 2.004602204236924;
else
    validateattributes(gValue,{'numeric'},{'scalar','nonnegative'},mfilename,'gValue')
end
validateattributes(DistanceAxis,{'numeric'},{'nonempty','increasing','nonnegative'},'getKernel','DistanceAxis')
validateattributes(TimeAxis,{'numeric'},{'nonempty','increasing','nonnegative'},'getKernel','TimeAxis')
validatestring(KernelBType,validKernelBTypes);
checklengths(TimeAxis,Background);

%--------------------------------------------------------------------------
%Kernel construction
%--------------------------------------------------------------------------
%Set digit precision to match accuracy of Fresnel integrals
OldPrecision = digits(15);

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

%If given, account for limited excitation bandwidth
if ~isempty(ExcitationBandwidth)
    wdd = exp(-wdd.^2/ExcitationBandwidth^2).*wdd;
end

Kernel = zeros(length(wdd));
%Compute dipolar kernel
for i=1:length(OvertoneCoeffs)
    %Allocate products for speed
    wddt = i*wdd.*TimeAxis;
    kappa = sqrt(6*wddt/pi);
    %Compute Fresnel integrals of 0th order
    C = fresnelC(kappa);
    S = fresnelS(kappa);
    Kernel = Kernel + OvertoneCoeffs(i)*sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S);
end
%Replace undefined Fresnel NaN value at time zero
Kernel(:,1) = 1;

%Build the background into the kernel
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
Kernel = Kernel.*repmat(Background,Dimension,1);

%Transpose
Kernel = Kernel';

%Return to old precision settings
digits(OldPrecision);
%--------------------------------------------------------------------------
return