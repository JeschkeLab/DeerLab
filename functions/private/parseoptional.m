function varargout = parseoptional(Properties,varargin)

if length(Properties)==1 && ~iscell(Properties)
    Properties = {Properties};
end
Properties = lower(Properties);

varargout = cell(length(Properties),1);

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
for i = 1:2:length(varargin(:))
    currentProperty = lower(varargin{i});
    if isa(currentProperty,'char')
        if any(strcmp(currentProperty,Properties))
            argoutidx = strcmp(currentProperty,Properties);
            varargout{argoutidx} = varargin{i+1};
        else
            warning('DA:parseoptional','''%s'' is not a valid property name.',currentProperty)
        end
    else
        error('The input is not a valid Name-Value pair.')
    end
    
end


end