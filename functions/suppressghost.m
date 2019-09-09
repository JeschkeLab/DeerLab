%
% SUPRESSGHOST Ghost distance suppression in multi-spin systems
%
%   cS = SUPRESSGHOST(S,n)
%   Suppresses multi-spin contributions to the signal (S) by means of the
%   power scaling approximation. The scaling is determined by the number of
%   radicals (n) in the system. The function returns the power scaled
%   signal (cS) without further normalization.
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function S = suppressghost(S,NRadicals)

if ~isreal(S)
    error('Input signal cannot be complex.')
end

if nargin<3
    NRadicals = 2;
end

Scaling = 1/(NRadicals - 1);

S = S.^Scaling;

end





