%
% MULTISTARTS Global optimization  - multiple starting point generator
%
%   [x0] = MULTISTARTS(n,x0,lb,ub)
%   Generates (n) starting points for non-linear constrained global
%   optimization. For each starting point a local minimum can be found. The
%   function uses a fixed random number generator to pick the starting
%   points from within the boundaries defined by the upper bounds (ub) and
%   lower bounds (lb).
%   The initial start point (x0) passed by the user is also included in the
%   set and returned with n-1 additional points.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function x0 = multistarts(n,x0,lb,ub)

if n<0
   error('The number of requested starting points must be n>0.') 
end

if numel(x0) ~= numel(lb) || numel(x0) ~= numel(ub)
   error('The lower/upper bound size(s) are not compatible with the initial guess vector x0.') 
end

% Ensure the use of row vectors
x0 = x0(:).';
lb = lb(:).';
ub = ub(:).';

% Generate n-1 new starting points within the bounds
x0 = [x0; (ub-lb).*rand(n-1,numel(x0)) + lb];

end






