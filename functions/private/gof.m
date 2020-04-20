%
% GOF Goodness of Fit statistics
%
%   stats = GOF(V,Vfit,Ndof,sigma)
%   Given some input data (V) and its corresponding fit (Vfit), with (Ndof)
%   degrees of freedom and with an underlying error with standard deviation
%   (sigma), returns a structure (stats) with multiple statistical
%   indicators of goodness of fit.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function stats = gof(V,Vfit,Ndof,sigma)

if nargin<4
    sigma = std(V - Vfit);
end
    
% Ensure correct dimensionality
V = V(:);
Vfit = Vfit(:);

% Get number of variables
N = numel(V);
% Extrapolate number of parameters
Q = Ndof - N;

% Reduced Chi-squared test
stats.chi2red = 1/Ndof*norm(V - Vfit)^2/sigma^2;

% R-squared test
stats.R2 = 1 - sum((V-Vfit).^2)/sum((Vfit-mean(Vfit)).^2);

% Root-mean square deviation
stats.RMSD = sqrt(sum((V-Vfit).^2)/N);

loglike = N*log(norm(V - Vfit)^2/N);
% Akaike information criterion
stats.AIC =  loglike + 2*Q;

% Corrected Akaike information criterion
stats.AICc = loglike + 2*Q + 2*Q*(Q+1)/(N-Q-1);

% Bayesian information criterion
stats.BIC =  loglike + Q*log(N);

end