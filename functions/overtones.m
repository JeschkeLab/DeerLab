%
% OVERTONES Analytical overtone coefficients of RIDME experiments
%
%   c = OVERTONES(n,Tmix,T1)
%   Computes the overtone/harmonics coefficients up to n-th order according to
%   analytical equations for a given mixing time Tmix and longitudinal 
%   relxation time T1. The function returns a n-element array containing
%   the coefficients.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

% see Keller et al, Phys.Chem.Chem.Phys., 2017, 19, 17856
% https://doi.org/10.1039/c7cp01524k
% (In that paper, overtone coefficients are indicated by P_k.)

function c = overtones(n,Tmix,T1)

% Validate input
if nargin<3 
   error('Not enough input arguments. Three are required.') 
end
validateattributes(n,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'n')
validateattributes(Tmix,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'Tmix')
validateattributes(T1,{'numeric'},{'scalar','nonempty','nonnegative'},mfilename,'T1')

% Determine overtone coefficients from the kinetic equations
c = 1 - exp(-Tmix/T1./(1:n));

% Normalize probability distribution
c = c/sum(c);

end
