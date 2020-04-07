%
% PARSEOPTIONAL Parser of property-value pair inputs
%
%   varargout = PARSEOPTIONAL(props,varargin)
%   Takes all inputs in (varargin) and identifies the allowed properties
%   specified by the cell array (props).
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function varargout = parseoptional(Properties,varargin)

%If only one property passed put it into a cell array
if length(Properties)==1 && ~iscell(Properties)
    Properties = {Properties};
end

% Ensure that all properties are lower case
Properties = lower(Properties);

% Check if flag is passed to ingore parsening errors
if any(contains(lower(Properties),'internal::parselater'))
    showErrors = false;
    Properties(contains(lower(Properties),'internal::parselater')) = [];
else
    showErrors = true;
end

% Pre-allocate output container
varargout = cell(length(Properties),1);

% If no varargins passed, just return
if isempty(varargin{1}), return; end

if length(varargin)==1
    varargin = varargin{1};
    if all(cellfun(@length,varargin)>1) && all(cellfun(@(x)isa(x,'cell'),varargin))
        varargin = varargin{1};
    end
end
if isempty(varargin{1})
    varargout = cell(nargout,1);
    return
end

% Parse property-value pairs in varargin
% --------------------------------------
for i = 1:2:length(varargin(:))
    currentProperty = lower(varargin{i});
    if isa(currentProperty,'char')
        
        % If property does exist then add the value to output list
        if any(strcmp(currentProperty,Properties))
            argoutidx = strcmp(currentProperty,Properties);
            varargout{argoutidx} = varargin{i+1};
            
        % If a property does not exist then show error unless supressed
        elseif showErrors
            % Get call stack
            parent = dbstack(1);
            % Compute Levenhstein distance between error property and rest
            levdist = cellfun(@(c)levendist(c,currentProperty),Properties);
            [mindist,idx] = min(levdist);
            % If the strings are not too different, suggest a correction
            if mindist<5 
            Correction = Properties{idx};
            error('DeerLab:parseoptional','There is no ''%s'' property on the %s function \n\nDid you mean: ''%s''',currentProperty,parent(1).file(1:end-2),Correction)
            else
            error('DeerLab:parseoptional','There is no ''%s'' property on the %s function',currentProperty,parent(1).file(1:end-2))
            end
        end
    else
        %If the varargin does not contain Property-value pairs, then show error
        error('DeerLab:parseoptional','The input is not a valid Property-Value pair.')
    end
    
end

end