%
% PERCENTILE Calculates the n-th percenitle of an array 
%            (similar to prctile function in Statistics Toolbox)
%
%   Y = PERCENTILE(X,p,dim)
%   Computes the p-th percentile (scalar value 0-100) of the
%   multi-dimensional array (X) along the dimension (dim).
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function Y = percentile(X,p,dim)

% Set requested dimension as the first dimension
dimIdx = 1:ndims(X);
dimIdx = dimIdx(dimIdx~=dim);
X = permute(X,[dim dimIdx]);

% Get size of data
sizeX = size(X);

% Vectorize all other dimensions
if numel(sizeX)>2
    X = reshape(X,[sizeX(1),prod(sizeX(2:end))]);
end

N = size(X,1);
% Sort data about first dimension
X = sort(X,1);
% Get list of available percentiles
pList = 100*(0.5:1:N-0.5)/N;
% Interpolate from list to requested percentile
Y = interp1(pList,X,p,'linear','extrap');

if numel(sizeX)>2
    % Reshape results back to original size
    Y = reshape(Y,sizeX(2:end));
end

end
