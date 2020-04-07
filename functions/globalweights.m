%
% GLOBALWEIGHTS Weights for global fitting of several signals
%
%   w = GLOBALWEIGHTS({S1,S2,S3,...})
%   Returns an N-element array of weights from the cell array (S) containing N
%   signals {S1,S2,S3,...,SN}. The weights are computed according to the
%   contribution of each signal to the ill-posedness of the inverse problem.
%
%   w = GLOBALWEIGHTS({S1,S2,S3,...},[n1,n2,n3,...])
%   The noise levels (n1,b2,n3,...,nN) for each signal can be passed as 
%   an N-element array.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.


function weights = globalweights(S,NoiseLevel)

% Get length of each signal
N = cellfun(@length,S);

% Estimate noise level of each signal if not provided
if nargin<2 || isempty(NoiseLevel)
    NoiseLevel = cellfun(@noiselevel,S);
end

if numel(S)~=numel(NoiseLevel)
  error('Number of signals and number of elements in noise level array must match.');
end

%Ensure use of column vectors
NoiseLevel = NoiseLevel(:);
N = N(:);

% Compute the weights
weights = 1./(N.*NoiseLevel);

% Normalize
weights = weights/sum(weights);

end
