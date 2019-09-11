%
% OVERTONES Analytical overtone coefficients of RIDME experiments
%
%   c = OVERTONES(n,Tm,T1)
%   Computes the overtone coefficients up to (n)-th order according to
%   analytical equations for a given mixing time (Tm) and longitudinal 
%   relxation time (T1). The function returns a n-element array containing
%   the coefficients.
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function OvertoneCoefficients = overtones(Order,MixingTime,T1)

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