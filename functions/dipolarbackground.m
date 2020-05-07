%
% DIPOLARBACKGROUND Multipathway background generator
%
%   B = DIPOLARBACKGROUND(t)
%   B = DIPOLARBACKGROUND(t,pathinfo,Bmodel)
%   B = DIPOLARBACKGROUND(___,'Property',value,___)
%
%  Computes the N-point multipathway background B from a N-point time axis t
%  (in microseconds), the dipolar pathways defined in pathinfo and a basis 
%  function Bmodel.
%  Optional arguments can be specified by name-value pairs.
%
%  Inputs:
%     t         N-element time axis, in microseconds
%     pathinfo  px2 or px3 array of modulation depths lambda, refocusing points
%               T0, and harmonics n for multiple pathways, each row contains
%               [lambda T0 n] or [lambda T0] for one pathway. If n is not given
%               it is assumed to be 1.
%     Bmodel    Function handle for a background model: @(t)bg_model(t,par)
%
%  Outputs:
%     B         N-element total multipathway background 
%
%  Name-value pairs:
%
%   'OvertoneCoeffs' - 1D array of coefficients for overtones in RIDME signals
%   'Renormalize'    - Re-normalization to ensure the equality B(0) == 1 is satisfied.
%                      (default = true)
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function B = dipolarbackground(t,pathinfo,Bmodel,varargin)

if nargin<3
   error('At least three inputs required: dipolarbackground(t,pathinfo,Bmodel)') 
end

% Check if user requested some options via name-value input
[OvertoneCoeffs,Renormalize] = ...
    parseoptional({'OvertoneCoeffs','Renormalize'},varargin);

if isempty(Renormalize)
    Renormalize = true;
else
    validateattributes(Renormalize,{'logical'},{'nonempty'},mfilename,'Renormalize');
end
if isempty(OvertoneCoeffs)
    OvertoneCoeffs = 1;
end
validateattributes(OvertoneCoeffs,{'numeric'},{'nonnegative'},mfilename,'OvertoneCoeffs');

if ~isa(Bmodel,'function_handle') || nargin(Bmodel)~=2 
    error(['For a model with multiple modulated pathways, B must be a ',...
           'function handle of the type: @(t,lambda) bg_model(t,par,lambda)']);
end

% Make sure all vectors are column vectors
t = t(:);
validateattributes(t,{'numeric'},{'nonempty'},mfilename,'t');

% Validation of the multi-pathway parameters
%-------------------------------------------------------------------------------
if ~isnumeric(pathinfo) || ~isreal(pathinfo)
    error('lambda/pathinfo must be a numeric array.');
end

if numel(pathinfo)==1
    lambda = pathinfo;
    pathinfo = [1-lambda NaN; lambda 0];
end

if ~any(size(pathinfo,2)==[2 3])
  error('pathinfo must be a numeric array with two or three columns.');
end
if any(isnan(pathinfo(:,1)))
  error('In pathinfo, NaN can only appear in the second column (refocusing time) e.g. path(1,:) = [Lam0 NaN];');
end

%Nomalize the pathway amplitudes to unity
pathinfo(:,1) = pathinfo(:,1)/sum(pathinfo(:,1));
lambda = pathinfo(:,1);
T0 = pathinfo(:,2);
if size(pathinfo,2)==2
  n = ones(size(T0));
else
  n = pathinfo(:,3);
end

% Combine all unmodulated components into Lambda0, and eliminate from list
unmodulated = isnan(T0);
Lambda0 = sum(lambda(unmodulated));
lambda(unmodulated) = [];
T0(unmodulated)  = [];
n(unmodulated) = [];

% Fold overtones into pathway list
nCoeffs = numel(OvertoneCoeffs);
if nCoeffs>0
  lambda = reshape(lambda*OvertoneCoeffs(:).',[],1);
  T0 = reshape(repmat(T0,1,nCoeffs),[],1);
  n = reshape(n*(1:nCoeffs),[],1);
end

nModPathways = numel(lambda);

% Construction of multi-pathway background function 
%-------------------------------------------------------------------------------
Bnorm = 1;
B = 1;
for p = 1:nModPathways
    B = B.*Bmodel(n(p)*(t-T0(p)),lambda(p));
    Bnorm = Bnorm.*Bmodel(-T0(p)*n(p),lambda(p));
end
if Renormalize
    B = B./Bnorm;
end

end