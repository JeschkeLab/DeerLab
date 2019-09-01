%
% GLOBALWEIGHTS Weights for global fitting of several signals
%
%   w = GLOBALWEIGHTS({S1,S2,S3,...})
%   Returns a N-point array of weights from the cell array (S) containing N
%   signals (S1,S2,S3,...,SN). The weights are computed according to the
%   contribution of each signal to the ill-posedness of the inverse problem.
%
%   w = GLOBALWEIGHTS({S1,S2,S3,...},[n1,n2,n3,...])
%   The noise levels (n1,b2,n3,...,nN) for each signal can be passed as 
%   a N-point array.
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function weights = globalweights(S,NoiseLevel)

%Get length of each signal
N = cellfun('length',S);
%If not passed by the user...
if nargin<2|| isempty(NoiseLevel)
    %...estimate the noise level for each signal
    NoiseLevel = cellfun(@noiselevel,S);
end
%Compute the weights
weights = sum(N.*NoiseLevel)./(N.*NoiseLevel);
%Normalize the weights
weights = weights/sum(weights);

end