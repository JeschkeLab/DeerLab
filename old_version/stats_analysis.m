function out = stats_analysis(handles,rms)
%STATS_ANALYSIS
%
% Written by H.C. Hyde, 2010 (Last updated: 2012-02-10)

% Set the one-tailed probability level P below (same as significance level ?)
% Use P<=0.32; default value: P=0.05 (~2 sigma)
%
% Approx values: P=0.32  (~1 sigma), P=0.05  (~2 sigma), P=0.01  (~2.6 sigma), P=0.003 (~3 sigma)
% Exact values : P=0.3173 (1 sigma), P=0.0455 (2 sigma), P=0.0124 (2.6 sigma), P=0.0027 (3 sigma)
% To calculate exact P-values: P = 1-erf(n/sqrt(2)), for any "n" sigma level

if exist('finv') ~= 2,
    out = [];
    return;
end;
  
P = 0.05;   %probability P / significance level ?

% Get number of fitted parameters and degrees of freedom for F statistic calculation
nobs = length(handles.A_cluster);   % # of data observations
if handles.model_updated
  switch handles.user_model
    case {'Two_Gaussians','Two_Rice3d'}
      nfp = length(find(handles.model_mask(1:5)));   % # of fitted parameters (exclude constraint indices)
    case 'Two_Gaussians_hom'
      nfp = length(find(handles.model_mask(1:6)));   % # of fitted parameters (exclude constraint indices)
    otherwise
      nfp = length(find(handles.model_mask));        % # of fitted parameters
  end
  ndf = nobs-nfp;   % # of degrees of freedom: sum of all data observation points minus # of fitted parameters
end

if handles.Tikh_updated
    if length(handles.Tikh_nfp_vector)>1
        %L-curve Tikhonov fit
        Iregpar = find(handles.rpv==handles.regpar);  % find best regularization parameter index
        nfp = handles.Tikh_nfp_vector(Iregpar);       % # of Tikhonov effective free parameters at best Reg. par.
    else
        %Single Tikhonov fit
        nfp = handles.Tikh_nfp_vector;   % # of Tikhonov effective free parameters at selected Reg. par.
    end
    ndf = nobs-nfp;   % # of Tikhonov effective degrees of freedom 
    %Note: From FTIKREG output, # of fitted parameters (nfp) equals the # of active constraints - 1)
end
 


% FOR COMPARISON OF MULTIPLE NON-NESTED (OR NESTED) MODELS
% --------------------------------------------------------
% Use the small-sample corrected Akaike Information Criterion (AICc) for 
% comparison of non-nested models fit to the same dataset. The AICc is used to
% compare the relative likelihood of 2 or more models, however it does not accept
% or reject a model at a statistical significance level and there is no reported 
% P value. The model with the lowest AICc value is most likely to be correct. 
% Given all AICc scores, it is simple to compute the probability (Akaike weight)
% that each of the tested models is correct. The AICc comparison is advantageous
% because there is no requirement for the models to be nested (which allows 
% comparison between a Tikhonov regularization and 2-gaussian model for example), 
% and more than 2 different models can be compared simultaneously. AICc is more 
% appropriate than the standard AIC when the sample size N is small or K is large
% (approx. N/K < 40), however AICc converges to AIC as N/K increases. Since most 
% fits of dipolar evolution data will be below this threshold, we only report the 
% more appropriate AICc value. 
%
% The standard AIC is given by the following general equation which includes the 
% maximum likelihood estimate L(theta_hat|data) and K = the number of estimatable
% parameters, which includes +1 for estimation of L(theta|data) at the optimal
% parameter set given by theta_hat.
%
% AIC = -2*log(L(theta_hat|data)) + 2*K
%
% For the least squares case, and assuming normally distributed errors with
% a constant variance, the AIC value is given by the following equation where 
% sigma^2_hat = sum(residuals^2)/N, and N is the number of independent observations:
%
% AIC = N*log(sigma^2_hat) + 2*K
%
% The above assumptions appear to be met by observation of fit residuals obtained 
% by fits of dipolar evolutions using the DeerAnalysis program. The readily available 
% goodness-of-fit output of the program is the root mean square (RMS) value, given by 
% RMS = sqrt(sum(residuals^2)/(N-1))). Therefore, sigma^2_hat ~ RMS^2 with the slight
% discrepancy of N-1 vs. N used for the RMS calculation in the DeerAnalysis program. 
% To be exact in the AIC calculation, we correct for the N-1 term in the RMS value 
% by using S = sigma^2_hat = RMS^2*(N-1)/N, therefore we have:
%
% AIC = N*log(S) + 2*K 
%     = N*log(RMS^2*(N-1)/N) + 2*K
%
% The small-sample corrected AICc value is calculated by adding a correction term:
%
% AICc = AIC + 2*K*(K+1)/(N-K-1)

N = nobs;    %number of observations
K = nfp+1;   %number of fitted parameters +1 to account for estimation of the RMS
S = rms^2*(N-1)/N;   %maximum likelihood estimate of sigma^2 for least squares fit
AIC  = N*log(S) + 2*K;            %AIC form for least squares case
AICc = AIC + 2*K*(K+1)/(N-K-1);   %corrected AIC (small-sample bias-adjusted AIC)

% Once multiple model fits and corresponding AICc values are obtained, for 
% comparative analysis we first examine the AICc differences: 
%
% delta_i = AICc_i - AICc_min, where 'i' is the model index and 'AICc_min' is the best model 
%
% The absolute AICc values are not informative, however the relative AICc differences are 
% very useful for model comparison and ranking, as well as for computing the positive
% "Akaike weights" and evidence ratios. The Akaike weights 'wi' give the relative likelihood
% of each model in the entire set of R models, which provides an intuitive probability-based 
% measure of relative likelihood. Note that the Akaike weights sum to 1 and are only valid
% for the set of R models and must be recomputed if the number of models changes. It also
% requires that one of the models must be the K-L best model of that set of R models, where
% K-L is the Kullback-Leibler information loss or distance (and whose minimization forms the
% basis of AIC theory).
%
% The likelihood of a given model 'gi' given data 'x' is:
%
% L(gi|x) ~ exp[-1/2*delta_i]
%
% The Akaike weights are normalized over the set of R models as:
%
% wi = exp[-1/2*delta_i] / sum(exp[-1/2*delta_r]), with summation limits r=1:R
%
% The evidence ratio is a measure of evidence for one model to be better than another
% in the K-L information sense and is given by the following equation. It can be used
% to compare any pair of models within the set of R models, but is often used to compare
% the best model (i=1) to other models j.
%
% L(gi|x)/L(gj|x) or equivalently (wi/wj)
%
% In conclusion, we reiterate the fact that the AICc comparisons presented above cannot
% be used for conventional hypothesis testing where a decision is to be reached, but rather 
% provide quantitative evidence for decision making in selection of the best model. The user
% should combine the quantitative evidence given by AIC theory with knowledge of the models
% and system under study to reach a conclusion as to which model is best (if any).
%
% A useful tabulation of all AIC results will look similar to this format:
%
% Model   delta_i   L(gi|x)    Akaike weight wi    Evidence Ratio w1/wi
% _____________________________________________________________________
%  1      0         1          0.6593               
%  2      1.4       0.4966     0.3274              2.0137
%  3      7.8       0.0202     0.0133              49.5714
%
% For reference:
% 1) Burnham, K.P. and Anderson, D.R. Model Selection and Multimodel
% Inference: A Practical Information-Theoretic Approach, 2nd Ed., Springer,
% New York (2002), 49-88.
%
%
%
% FOR COMPARISON OF CONSTRAINED VS. UNCONSTRAINED FIT WITH SAME MODEL
% -------------------------------------------------------------------
% Use the likelihood-ratio method as a one-tailed F-test to determine if the 
% constrained model parameter set is equal to the unconstrained parameter
% set at probability level P. The following equation describes the total
% variance of the fit in terms of RMS error:
%
% RMS(thresh) = RMS(min)*sqrt(1+p/(N-p)*F(p,N-p,1-?)) 
%
% where F is the upper 1-? quantile of the central Fisher’s F distribution 
% evaluated at p and N-p degrees of freedom; p is the number of simultaneously
% fitted parameters; N is the number of fitted data points. See [Seber&Wild 1989]
% for reference. The significance level ? is typically set at ?=0.05 (2?). The 
% probability level P, same as ?, is user-defined at the top of this m-file.
%
% F(p,N-p,1-?) is evaluated in Matlab using the inverse of the F cumulative
% distribution function as follows: F(p,N-p,1-?) = finv(1-?,p,N-p).
%
% RMS(thresh) represents a p-dimensional surface contour that is the boundary of 
% the 100(1-P)% confidence region. By the duality of confidence intervals 
% and hypothesis tests, we construct an F-test under the following 
% hypothesis: Ho: ?=?o vs. Ha: ???o, where ? and ?o are the parameter sets
% obtained by constrained and unconstrained model fits, respectively.
% 
% The F-test: If RMS(constrained) <= RMS(thresh), we accept the null hypothesis 
% that the unconstrained and constrained models have equivalent parameter sets at
% the probability level P. If RMS(constrained) > RMS(thresh), we accept the 
% alternative hypothesis that the unconstrained and constrained models have 
% significantly different parameter sets. If the unconstrained model fits the 
% dipolar evolution well, but the constrained model reaches significance, this
% suggests that the nonlinear constraint is not valid for this dataset, 
% this model, etc.
%
% Confidence intervals: The 100(1-?)% confidence interval of a selected parameter
% is obtained by projection of the 100(1-?)% confidence region onto the selected
% parameter axis. If a confidence region has not been calculated, the confidence 
% interval of a selected parameter can be calculated by finding its lower and upper
% values that satisfy the likelihood-ratio, e.g. RMS(thresh), at the desired 
% significance level. This method involves examination of the error surface, 
% and is also referred to as a support plane analysis [Straume&Johnson 1991].
% The procedure is to first obtain the optimal parameter set from an unconstrained
% fit with minimum error, e.g. RMSD(min). Next, the selected parameter is incremented
% away from its optimal value and the least-squares fit is repeated with the parameter
% fixed at the new value. This process continues until the fit error exceeds the error
% threshold value given by the likelihood-ratio. The parameter value that defines
% the exact confidence interval boundary is then obtained by interpolation (cubic) 
% of several test parameter values near the error threshold. This procedure is then
% repeated in the opposite direction to obtain the opposite confidence interval bound.
% This same procedure can be used to construct the error surface for the nonlinear 
% parameter <r2>/<r1>. Simply fix the mean distance ratio k to selected values, run 
% the fit and record the RMS value. Repeat this procedure to build a set of values for
% RMS vs. k, which define the error surface for the mean distance ratio. 
% 
% For reference:
% 1) Seber, G.A.F. and Wild, C.J. (1989). Nonlinear Regression, John Wiley and Sons:
%    New York, (Chapter 5), pp. 191-269.
% 2) Straume, M.; Frasier-Cadoret, S.G.; Johnson, M.L. Least-squares analysis of fluorescence
%    data. In Topics in fluorescence spectroscopy, Vol. 2: Principles, Ed JR Lakowicz. 
%    Plenum Press, New York 1991, 177-240.
%
% Note that the constrained model has nonlinear parameter constraints and the 
% same number of fitted parameters. This F-test is only valid for comparison of an
% unconstrained (best-fit) and constrained fit of the SAME model and SAME dataset.
Fndf_numer = nfp;   % # of degrees of freedom in F distribution numerator
Fndf_denom = ndf;   % # of degrees of freedom in F distribution denominator
Fpct = finv(1-P,Fndf_numer,Fndf_denom);   %F percentile for comparison of rms values 


% Define output results
out.rms        = rms;          % RMS value of current fit 
out.rmsT       = sprintf('%9.7f',out.rms);   %text string of RMS
out.nobs       = nobs;         % # of data observations
out.nfp        = nfp;          % # of fitted parameters
out.ndf        = ndf;          % # of degrees of freedom
out.Fndf_numer = Fndf_numer;   % F distribution degrees of freedom in numerator
out.Fndf_denom = Fndf_denom;   % F distribution degrees of freedom in denominator
out.P          = P;            % statistical significance probability level (same as alpha)
out.alpha      = P;            % statistical significance level alpha (same as P)

% Calculate F percentile and RMSD threshold for model fits only.
%  Note: Does not apply to Tikhonov fits.
if handles.model_updated
    out.Fpct  = Fpct;   %F percentile for comparison of RMS values 
    out.FpctT = sprintf('%9.7f',out.Fpct);  %text string of F percentile
    out.rms_thresh  = rms*sqrt(1+nfp/ndf*Fpct);   %RMS statistical threshold value
    out.rms_threshT = sprintf('%9.7f',out.rms_thresh);  %text string of RMS threshold
end

% Calculate AICc for model fit or Tikhonov regularization
if handles.model_updated || handles.Tikh_updated
    out.AICc  = AICc;   %AICc = corrected AIC
    out.AICcT = sprintf('%9.7f',out.AICc);   %text string of AICc
end


