%
% NOISELEVEL Estimate the noise standard deviation in a signal
%
%   sigma = NOISELEVEL(V2D)
%   sigma = NOISELEVEL(Vco)
%   sigma = NOISELEVEL(V)
%   sigma = NOISELEVEL(V,filter)
%   sigma = NOISELEVEL(V,Vref)
%
%   Returns the standard deviation estimation of the noise in a given signal.
%
%   (V2D) is a 2D-dataset of different scans, the noise standard deviation
%   is estimated from the deviations between scans. The second dimension of
%   (V2D) must contain the different scans. The function returns the standard
%   deviation of the averaged signal not of the individual scans.
%
%   If a 1D signal (V) is given, the noise level is estimated via filtering
%   of the signal with a moving mean filter. The nature of the filter can
%   be specified by means of a string (filter).
%
%   If a reference model signal (Vref) is given, the noise level is
%   estimated from the difference between both signals.
%
%   If the input signal (Vco) contains an imaginary component, the noise
%   level is estimated form the imaginary component after phase optimization.
%
%   Inputs:
%       V2D         NxM-element datasets of single scans of a dipolar signal
%       Vco         N-element complex dipolar signal
%       V           N-element dipolar signal
%       filter      Type of filtering:
%                          - 'movmean'  Moving mean filter
%                          - 'savgol'   Savitzky-Golay filter
%       Vref        N-element reference signal
%
%   Outputs:
%       sigma       Estimated noise standard deviation

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function sigma = noiselevel(V,argin)


validateattributes(V,{'numeric'},{'nonempty'},mfilename,'V')

% Input: noiselevel(V2D)
if ~isvector(V)
    estimationMethod = '2D';
    if nargin>1
        error('For 2D-datasets, only one input is required.')
    end

% Input: noiselevel(V)
elseif isvector(V) && nargin<2 && isreal(V)
    V = V(:);
    estimationMethod = 'filtering';
    filterType = 'movmean';

% Input: noiselevel(Vco)
elseif isvector(V) && nargin<2 && ~isreal(V)
    V = V(:);
    estimationMethod = 'complex';
    
% Input: noiselevel(V,filter)
elseif isvector(V) && ischar(argin) && isreal(V)
    V = V(:);
    estimationMethod = 'filtering';
    filterType = argin;

% Input: noiselevel(V,Vref)
elseif isvector(V) && ~ischar(argin) && isreal(V)
    estimationMethod = 'reference';
    Vref = argin;
    V = V(:);
    Vref = Vref(:);
    validateattributes(Vref,{'numeric'},{'nonempty'},mfilename,'Vref')
    if numel(V) ~= numel(Vref)
       error('The input and reference signal must have the same number of elements.') 
    end
else
        error('The input is not valid.')
end

switch estimationMethod
    
    case '2D'
        % Estimate standard deviations for all time point, and average over scans
        if size(V,2)<10
            warning('Only a few scans are given. Noise standard deviation estimate will be inaccurate.');
        end
        sigma = std(V,0,2);
        sigma = mean(sigma);
        
    case 'filtering'
        
        switch filterType
            case 'movmean'
                Vfilt = movmean(V,12);
            case 'savgol'
                Vfilt = savGol(V,11,11,3);
        end
        sigma = std(V - Vfilt);
        
    case 'complex'
        [~,Vim] = correctphase(V);
        sigma = std(Vim);
        
    case 'reference'
        sigma = std(V - Vref);
end



end
