%
% NOISELEVEL Estimate the noise level (standard deviation)
%
%   STD = NOISELEVEL(S)
%   Returns the standard deviation estimation of the noise in the signal S.
%   The estimation is done from the last quarter of the signal. 
%
%   STD = NOISELEVEL(S,M)
%   Returns the standard deviation estimation of the noise in the signal S.
%   The estimation is done from the last M-points of the N-point signal. 
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function Level = noiselevel(Signal,M)

%Get signal length
N = length(Signal);
%Extract the piece of signal to estimate
if nargin<2 || isempty(M)
    Cutoff = Signal(ceil(4/5*N):N);
else
    Cutoff = Signal(ceil(N-M):N);   
end
%Estimate the noise level
Cutoff = Cutoff - mean(Cutoff);
Level = std(Cutoff);

end