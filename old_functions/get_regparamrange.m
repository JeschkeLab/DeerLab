% get_alpharange  Determine regularization parameter range
%
%  alpha = get_regparamrange(K,L)
%  alpha = get_regparamrange(K,L,noise)
%
%  Determines an appropriate range of possible values for the regulariation
%  parameter alpha, based on the GSVD of K and L.
%
%  Input:
%    K          kernel matrix, N x N
%    L          regularization matrix, M x N
%    noise_std  noise standard deviation of form factor (0 - 1) (default 0)
%
%  Output
%    alpha  vector of alpha values, in decreasing order

function alpha = get_regparamrange(K,L,noise_std)

if nargin<3
  noise_std = 0;
end

lgregpar_inc = 0.1;  % resolution of log10(alpha)

% Set alpha range
%-------------------------------------------------------------
minmax_ratio = 16*eps*1e6;  % ratio of smallest to largest alpha

% Scaling by noise. This improves L curve corner detection for DEER
minmax_ratio = minmax_ratio*2^(noise_std/0.0025);

% Get generalized singular values of K and L
%-------------------------------------------------------------
singval = gsvd(K,L,0);
singval = singval(1:end-2); % remove two inf (for L = second deriv)
singval = singval(end:-1:1); % sort in decreasing order

% Calculate range based on singular values
%-------------------------------------------------------------
lgalpha_max = log10(singval(1));
lgalpha_min = log10(max([singval(end),singval(1)*minmax_ratio]));
lgalpha_max = floor(lgalpha_max/lgregpar_inc)*lgregpar_inc;
lgalpha_min = ceil(lgalpha_min/lgregpar_inc)*lgregpar_inc;
lgalpha = lgalpha_max:-lgregpar_inc:lgalpha_min;
alpha = 10.^lgalpha;
