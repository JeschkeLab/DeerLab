% NORM_INV Inverse of the normal cumulative distribution function (cdf).
%
%   X = NORM_INV(p,mu,sigma)
%   Returns the inverse cdf for the normal distribution with mean (mu)
%   and standard deviation (sigma), evaluated at the values in (p).
%
%   Default values for (mu) and (sigma) are 0 and 1, respectively.
%
% Based on the comments of Wayne King on MATLAB Answers:
% https://ch.mathworks.com/matlabcentral/answers/82596-workaround-for-the-norminv-function-in-statistics-toolbox#answer_92266%
%
% This is a workaround for the norminv() function from the Statistics and
% Machine Learning Toolbox from MATLAB.

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.
function [z] = norm_inv(p,mu,sigma)

if nargin<2
    mu = 0;
    sigma = 1;
elseif nargin<3
    mu = 0;
end

if numel(mu)+numel(sigma)+numel(p)~=3
    error('Input arguments must be scalar values.')
end

z = sigma*(-sqrt(2)*erfcinv(2*p)) + mu;

end