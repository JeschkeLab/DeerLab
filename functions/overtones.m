%
% OVERTONES Analytical overtone coefficients of RIDME experiments
%
%   c = OVERTONES(n,Tm,T1)
%   Computes the overtone coefficients up to (n)-th order according to
%   analytical equations for a given mixing time (Tm) and longitudinal 
%   relxation time (T1). The function returns a n-element array containing
%   the coefficients.
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

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