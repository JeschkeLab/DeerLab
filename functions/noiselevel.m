%
% NOISELEVEL Estimate the noise level (standard deviation)
%
%   level = NOISELEVEL(S)
%   Returns the standard deviation estimation of the noise in the signal S.
%   The estimation is done from the last quarter of the signal.
%
%   level = NOISELEVEL(S,M)
%   Returns the standard deviation estimation of the noise in the signal S.
%   The estimation is done from the last M-points of the N-point signal.
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function Level = noiselevel(S,M)

validateattributes(S,{'numeric'},{'2d','nonempty'})
%Get signal length
N = length(S);

%Validate input
if nargin<2 || isempty(M)
    %Extract the piece of signal to estimate
    M = 1/5*N;
else
    validateattributes(M,{'numeric'},{'scalar','nonempty','nonnegative'})
end
if M>N
    error('Second argument cannot be longer than the length of the signal.')
end

%Extract the piece of signal to estimate
Cutoff = S(ceil(N-M):N);

%Estimate the noise level
Cutoff = Cutoff - mean(Cutoff);
Level = std(Cutoff);

end