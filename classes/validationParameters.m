classdef validationParameters
    
    
    properties
        name
        trials
        vector
        range
        sampling = 'lin';
    end
    
    methods
        
        function obj = set.name(obj,string)
            validateattributes(string,{'char'},{})
            obj.name = string;
        end
        
        function obj = set.range(obj,array)
            validateattributes(array,{'numeric'},{'2d'})
            if length(array)~=2
               error('Property ''range'' must be an array of two elements.') 
            end
            obj.range = array;
        end
        
        function obj = set.trials(obj,value)
            validateattributes(value,{'numeric'},{'scalar','nonnegtive'})
            obj.trials = value;
        end
        
        function obj = set.vector(obj,array)
            validateattributes(array,{'numeric','cell'},{'2d'})
            obj.vector = array;
            if isempty(obj.trials)
                obj.trials = length(array);
            end
        end
        
        function obj = set.sampling(obj,string)
            validateattributes(string,{'char'},{'nonempty'})
            validstrings = {'lin','rand'};
            validatestring(string,validstrings);
            obj.sampling = string;
        end
        
    end
    
    
    
    
    
    
    
end