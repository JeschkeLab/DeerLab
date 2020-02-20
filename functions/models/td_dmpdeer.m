%
% TD_DMPDEER Time-domain dipolar multi-pathway (DMP) model for multi-pulse DEER
%
%     V = td_dmpdeer(t,r,P,pathways,B)
%
% Inputs:
%   t         time axis (us)
%   r         distance axis (nm)
%   P         distance distribution (nm^-1)
%   pathinfo  amplitudes and refocusing times of dipolar pathways, one pathway per row
%   B         background signal; or background model function if more than one
%             pathway are given
%
% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function V = td_dmpdeer(t,r,P,pathinfo,B)

if nargin<4
   error('Not enough input arguments.') 
end
if nargin<5
  B = [];
end

% Input validation
validateattributes(t,{'numeric'},{'nonempty'},mfilename,'t');
validateattributes(r,{'numeric'},{'nonnegative','nonempty'},mfilename,'r');
validateattributes(P,{'numeric'},{'nonnegative','nonempty'},mfilename,'P');
validateattributes(pathinfo,{'numeric'},{'nonempty'},mfilename,'pathinfo');
if ~isempty(B)
    validateattributes(B,{'numeric','function_handle'},{},mfilename,'B');
    if isnumeric(B) 
        checklengths(t,B);
    end
end

K = dipolarkernel(t,r,pathinfo,B);
V = K*P(:);
