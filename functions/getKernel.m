function [Kernel] = getKernel(TimeAxis,DistanceAxis,Background,varargin)

%--------------------------------------------------------------------------
%Input parsing
%--------------------------------------------------------------------------
if nargin<3 || isempty(Background)
   Background = ones(1,length(TimeAxis));  
end
validKernelBTypes = ["none","full","sqrt"];
[KernelBType,ExcitationBandwidth] = parseOptional({'KernelBType','ExcitationBandwidth'},varargin);
if isempty(KernelBType)
    KernelBType = 'sqrt';
end
if ~isempty(ExcitationBandwidth)
    validateattributes(ExcitationBandwidth,{'numeric'},{'scalar','nonnegative'})
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
validateattributes(DistanceAxis,{'numeric'},{'nonempty','increasing','nonnegative'},'DistanceAxis')
validateattributes(TimeAxis,{'numeric'},{'nonempty','increasing','nonnegative'},'TimeAxis')
validatestring(KernelBType,validKernelBTypes);
checklengths(TimeAxis,Background);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%Kernel construction
%--------------------------------------------------------------------------
%Set digit precision to match accuracy of Fresnel integrals
OldPrecision = digits(15);

Dimension = length(TimeAxis);
TimeStep = mean(diff(TimeAxis));

%Numerical dipolar frequency at 1 nm for g=ge
ny0 = 52.04;
%Convert time step to miliseconds if given in nanosecond
if TimeStep>1
    TimeStep = round(TimeStep)/1000;
    TimeAxis = round(TimeAxis)/1000;
end

%Numerical angular dipolar frequency at 1 nm for g=ge
w0 = 2*pi*ny0;
%Get vector of dipolar frequencies
wdd=w0./(DistanceAxis.^3);
wdd = wdd;

%If given, account for limited excitation bandwidth
if ~isempty(ExcitationBandwidth)
    wdd = exp(-wdd.^2/ExcitationBandwidth^2).*wdd;
end

%Allocate products for speed
wddt = wdd.*TimeAxis;
kappa = sqrt(6*wddt/pi);
%Compute Fresnel integrals of 0th order
C = fresnelC(kappa);
S = fresnelS(kappa);

%Compute dipolar kernel
Kernel = sqrt(pi./(wddt*6)).*(cos(wddt).*C + sin(wddt).*S);
%Replace undefined Fresnel NaN value at time zero
Kernel(:,1) = 1;

%Build the background into the kernel
switch KernelBType
    case 'none'
        Background = ones(Dimension,1);
    case 'sqrt'
        Background = sqrt(Background);
end
Kernel = Kernel + (1/Background(1) - 1);
Kernel = Kernel.*repmat(Background,Dimension,1);


%Transpose
Kernel = Kernel';

%Return to old precision settings
digits(OldPrecision);
%--------------------------------------------------------------------------
return