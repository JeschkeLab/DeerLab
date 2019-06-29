classdef daopts
    
    properties
        BackgroundModel = 'exponential'
        FitStartSearch = 'auto'
        FitStartPos = [];
    end
    
    properties (SetAccess = private)
        
    end
    
    properties (Hidden)
        ID
    end
    
    methods
        
        function obj = daopts(varargin)
            if nargin>0
                if mod(nargin,2)
                    error('Error \n Wrong number of arguments')
                end
                for i=1:2:length(varargin)
                    obj = setProperty(obj,varargin{i},varargin{i+1});
                end
            end
            JavaID =  java.util.UUID.randomUUID;
            obj.ID = JavaID.toString;
        end
        
        function obj = set.BackgroundModel(obj,string)
            allowedInput = {'exponential','polynomial','polyexp','fractal'};
            if any(strcmp(allowedInput,string)) && ~isa(string,'numerical')
                obj.BackgroundModel = string;
            else
                error('daopts:incorrectType','Error \n ''%s'' is not a valid input of the BackgroundModel property.',string)
            end
        end
        
        function obj = set.FitStartSearch(obj,string)
            allowedInput = {'auto','manual'};
            if any(strcmp(allowedInput,string)) && ~isa(string,'numerical')
                obj.BackgroundModel = string;
                if strcmp(string,'auto')
                    obj.FitStartPos = [];
                else
                    obj.FitStartPos = 1;
                end
            else
                error('daopts:incorrectType','Error \n ''%s'' is not a valid input of the BackgroundModel property.',string)
            end
        end
        
    end
    
    methods(Access = private)
        function obj = setProperty(obj,Property,Value)
            if ~isa(Property,'char') && ~isa(Property,'string')
                error('daopts:incorrectType','Error \n Invalid first data argument')
            end
            eval(sprintf('obj.%s = Value;',Property));
        end
    end
    
end