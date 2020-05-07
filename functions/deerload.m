%
% DEERLOAD  Load experimental EPR data
%
%   y = DEERLOAD(FileName)
%   [x,y] = DEERLOAD(FileName)
%   [x,y,Pars] = DEERLOAD(FileName)
%   [x,y,Pars,FileN] = DEERLOAD(FileName)
%   ... = DEERLOAD(FileName,Scaling)
%   ... = DEERLOAD
%
%   Read spectral data from a file specified in the string
%   'FileName' into the arrays x (abscissa) and y (ordinate) 
%   in microseconds. The structure Pars contains entries from the
%   parameter file, if present.
%
%   All strings in the parameter structure containing numbers
%   are converted to numbers for easier use.
%
%   If FileName is a directory, a file browser is
%   displayed. If FileName is omitted, the current
%   directory is used as default. DEERLOAD returns the
%   name of the loaded file (including its path) as
%   fourth parameter FileN.
%
%   Supported formats are identified via the extension
%   in 'FileName'. Extensions:
%
%     Bruker BES3T:        .DTA, .DSC
%     Bruker ESP, WinEPR:  .spc, .par
%     SpecMan:             .d01, .exp
%     Magnettech:          .spe (binary), .xml (xml)
%     Active Spectrum:     .ESR
%     Adani:               .dat, .json
%     JEOL:                (no extension)
%
%     MAGRES:              .PLT
%     qese, tryscore:      .eco
%     Varian:              .spk, .ref
%     ESE:                 .d00, .exp
%
%     For reading general ASCII formats, use textscan(...)
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function varargout = deerload(FileName,Scaling)

if nargout>4
    error('Please provide 1, 2, 3 or 4 output arguments!');
end

if nargin<1
    FileName = pwd;
end

if nargin<2
    Scaling = '';
end

% Prepare output container
varargout = cell(1,nargout);

% Load the data using eprload
[varargout{:}] = eprload(FileName,Scaling);

% Convert time axis to microseconds if loaded in nanoseconds
if length(varargout)>1
    t = varargout{1};
    V = varargout{2};
    if iscell(t)
        t = t{1};
        tmp{1} = V(:,1);
        tmp{2} = V(:,2);
        V = tmp;
    end
    
    t = t/1e3; % nanoseconds -> microseconds
    
    varargout{1} = t;
    varargout{2} = V;
end

return
