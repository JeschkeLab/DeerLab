% T_INV Inverse of the normal cumulative distribution function (cdf).
%
%   X = T_INV(p,v)
%   Inverse of Student's T cumulative distribution function (cdf).
%   Returns the inverse of Student's T cdf with (v) degrees
%   of freedom, at the values in (p)
%
% Based on the comments of Star Strider on MATLAB Answers:
% https://ch.mathworks.com/matlabcentral/answers/260564-how-to-implement-tinv-manually-without-statistics-and-machine-learning-toolbox#answer_203439
%
% This is a workaround for the norminv() function from the Statistics and
% Machine Learning Toolbox from MATLAB.

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function z = t_inv(p,v)

if numel(v)+numel(p)~=2
    error('Input arguments must be scalar values.')
end

tdist1T = @(t,v) 1-(1-(1-betainc(v/(v+t^2),v/2,0.5)))/2;
z = fzero(@(tval) (max(p,(1-p)) - tdist1T(tval,v)), 0);

%Determine sign of critical value
if p>0.5
    s = -1;
else
    s = +1;
end
z = s*z;

end
