%
% SUPRESSGHOST Ghost distance suppression in multi-spin systems
%
%   Vs = SUPRESSGHOST(V,n)
%   Suppresses multi-spin contributions to the signal (V) by means of the
%   power scaling approximation. The scaling is determined by the number of
%   radicals (n) in the system. The function returns the power scaled
%   signal (Vs).
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

function Vs = suppressghost(V,nSpins)

if nargin~=2
    error('Two inputs are required: the signal V, and the number of spins n.');
end

if ~isreal(V)
    error('Input signal cannot be complex.')
end

if numel(nSpins)~=1 || ~isreal(nSpins) || mod(nSpins,1)~=0
   error('Number of spins (second input) must be a positive integer.');
end

Scaling = 1/(nSpins-1);

Vs = V.^Scaling;

end
