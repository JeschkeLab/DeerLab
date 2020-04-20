%
% NOISELEVEL Estimate the noise standard deviation
%
%   level = NOISELEVEL(S)
%   Returns the standard deviation estimation of the noise in the signal S.
%   The estimation is done from the last fifth of the signal.
%
%   level = NOISELEVEL(S,M)
%   Returns the standard deviation estimation of the noise in the signal S.
%   The estimation is done from the last M points of the N-point signal.
%
%   level = NOISELEVEL(S)
%   If S is a 2D-dataset of different scans, the noise standard deviation 
%   is estimated from the deviations between scans. The second dimension of
%   V must contain the different scans. The function returns the standard 
%   deviation of the averaged signal not of the individual scans.  

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function level = noiselevel(V,M)

validateattributes(V,{'numeric'},{'2d','nonempty'})

if all(size(V)>1)
    is2D = true;
else
    is2D = false;
end


if is2D
    N = size(V,2);
    % Estimate the standard deviation of the single scans
    level = mean(std(V,[],2));
    % Estimate the standard deviation after averaging
    level = level/sqrt(N);
else
    V = V(:);
    
    % Get signal length
    N = length(V);
    
    % Validate input
    if nargin<2 || isempty(M)
        M = ceil(1/5*N);
    else
        validateattributes(M,{'numeric'},{'scalar','nonempty','nonnegative'})
    end
    if M>N
        error('Second argument cannot be longer than the length of the signal.')
    end
    
    % Extract the piece of signal to estimate
    idx = ceil(N-M):N;
    idx = idx(:);
    Cutoff = V(idx);
    
    % Fit a line to remove possible oscillation fragment
    opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);
    lineparam = fminsearch(@(x)norm(x(1) + x(2)*idx  - Cutoff)^2,rand(2,1),opts);
    linearfit = lineparam(1)+ lineparam(2)*idx;
    Cutoff = Cutoff - linearfit;
    
    % Estimate the noise level
    Cutoff = Cutoff - mean(Cutoff);
    level = std(Cutoff);
end


end
