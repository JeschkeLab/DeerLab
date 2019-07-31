function varargout = parseoptional(Properties,varargin)
if length(Properties)==1 && ~iscell(Properties)
    Properties = {Properties};
end



varargout = cell(length(Properties),1);

if ~isempty(varargin{1})
    if length(varargin)==1
        varargin = varargin{1};
    end
    for i = 1:2:length(varargin(:))
        currentProperty = varargin{i};
        if isa(currentProperty,'char')
            if any(strcmp(currentProperty,Properties))
            argoutidx = strcmp(currentProperty,Properties);
            varargout{argoutidx} = varargin{i+1};
            else
                error('''%s'' is not a valid property name.',currentProperty)
            end
        else
            error('The input is not a valid Name-Value pair.')
        end
        
    end
    
end

end