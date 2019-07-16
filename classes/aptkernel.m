classdef aptkernel < DeerAnalysis_object
    
    properties
        Base
        NormalizationFactor
        FreqAxis
        TimeAxis
        Crosstalk
    end
   
    methods
        
        function obj = aptkernel(varargin)
            if nargin>0
                if mod(nargin,2)
                    error('Wrong number of arguments')
                end
                for i=1:2:length(varargin)
                    obj = setProperty(obj,varargin{i},varargin{i+1});
                end
            end
        end
        
        function [Base,NormalizationFactor,FreqAxis,TimeAxis,Crosstalk] = dismountAPTkernel(obj)
            Base = obj.Base;
            NormalizationFactor = obj.NormalizationFactor;
            FreqAxis = obj.FreqAxis;
            TimeAxis = obj.TimeAxis;
            Crosstalk = obj.Crosstalk;
        end
    end
    
end