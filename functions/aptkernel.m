% 
% APTKERNEL Computes the dipolar interaction kernel and elements required 
%           for the approximate Pake transformation (APT).  
%
%       K = APTKERNEL(t) 
%       Computes a structure (K) containing the (N/2-2)xN point kernel,
%       the (N/2-2) point array of normalization factors, (N/2-2) point
%       frequency axis and the (N/2-2)x(N/2-2) crosstalk matrix corresponding
%       to the N-point time axis (t).
%
%       K = APTKERNEL(T,'ExcitationBandwidth',w)
%       The excitation bandwidth (w) of the experiment can be passed as an
%       option to account for it in the kernel.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

% Reference:
%    G. Jeschke, A. Koch, U. Jonas, A. Godt
%    Direct Conversion of EPR Dipolar Time Evolution Data to Distance Distributions
%    J.Magn.Reson. 155, 72-82 (2002)
%    https://doi.org/10.1006/jmre.2001.2498

function APTkernel = aptkernel(t,varargin)

% Check if user requested some options via name-value input
[ExcitationBandwidth] = parseoptional({'ExcitationBandwidth'},varargin);
if ~isempty(ExcitationBandwidth)
    validateattributes(ExcitationBandwidth,{'numeric'},{'scalar','nonnegative'})
end
if iscolumn(t)
   t = t.'; 
end
validateattributes(t,{'numeric'},{'nonempty','increasing'},'t')

% Use absolute time scale, required for negative times
t = abs(t);

% Turn off warnings to avoid ill-conditioned warnings 
warning('off','all')

FreqElement = 1/(2*max(t));
nFreq = floor(length(t)/2)-2;
FreqAxis = FreqElement*((1:nFreq)+1/4*ones(1,nFreq));

% Compute dipolar kernel using Fresnel integrals
wdd = 2*pi*FreqAxis(:);
ph = wdd.*t;
kappa = sqrt(6*ph/pi);
Base = (cos(ph).*fresnelC(kappa) + sin(ph).*fresnelS(kappa))./kappa;
Base(isnan(Base)) = 1; 

% If given, account for limited excitation bandwidth
if ~isempty(ExcitationBandwidth)
    Base = exp(-wdd'.^2/ExcitationBandwidth^2)'.*Base;
end

Base = Base*mean(diff(FreqAxis));

% Calculate normalization factors and cross talk matrix (eqns 19 and 20)
a = (Base.*t)*Base.';
NormalizationFactor = diag(a);
crosstalk = a./NormalizationFactor;

% Construct the kernel object to be passed later to the APT.m function
APTkernel = struct('Base',Base,...
                      'NormalizationFactor',NormalizationFactor,...
                      'FreqAxis',FreqAxis(:),...
                      't',t(:),...
                      'Crosstalk',crosstalk);

% Turn warnings back on
warning('on','all')

end
