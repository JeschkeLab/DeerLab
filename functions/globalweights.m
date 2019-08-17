%
% GLOBALWEIGHTS Weights for global fitting of several signals
%
%   w = GLOBALWEIGHTS({S1,S2,S3,...})
%   Returns a N-point array of weights from the cell array (S) containing N
%   signals (S1,S2,S3,...,SN). The weights are computed according to the
%   contribution of each signal to the ill-posedness of the inverse problem.
%
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function weights = globalweights(Signal)

N = cellfun('length',Signal);
NoiseLevel = cellfun(@noiselevel,Signal);
weights = sum(N.*NoiseLevel)./(N.*NoiseLevel);
weights = weights/sum(weights);

end