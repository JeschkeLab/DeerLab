function [nonlinearInequality,nonlinearEquality] = unitIntConstraint(x)
%-------------------------------------------
% Non-linear conditions for fmincon solver
%-------------------------------------------
%   /
%  |  nonlinearInequality(x) <= 0 for all x
% <
%  |  nonlinearEquality(x) == 0 for all x
%   \


%Nonlinear inequality not required so set to automatic fullfillment
nonlinearInequality = 0;

%Enforce unity integral of fitted distance distributions
nonlinearEquality = sum(x) - 1;

end