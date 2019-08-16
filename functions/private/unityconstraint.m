function [nonlinearInequality,nonlinearEquality] = unityconstraint(x,dr)
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
nonlinearEquality = sum(x)*dr - 1;

end