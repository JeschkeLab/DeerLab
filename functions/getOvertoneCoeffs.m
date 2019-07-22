function OvertoneCoefficients = getOvertoneCoeffs(Order,MixingTime,T1)

%Validate input
if nargin<3 
   error('Not enough input arguments') 
end
validateattributes(Order,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'Order')
validateattributes(MixingTime,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'MixinTime')
validateattributes(T1,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'T1')

%Determine overtone coefficients from the kinetic equations
OvertoneCoefficients = zeros(1,Order);
for i = 1:Order
   OvertoneCoefficients(i) = 1 - exp(-MixingTime/(i*T1));
end

%Normalize probability distribution
OvertoneCoefficients = OvertoneCoefficients/sum(OvertoneCoefficients);

end