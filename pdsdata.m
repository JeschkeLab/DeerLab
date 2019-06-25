classdef pdsdata
  
  properties
    TimeAxis
    RawData
    FormFactor
    DipEvoFctn
    ModDepth
  end
  
  properties (SetAccess = private)
    TimeStep
    Length
    TimeUnits
  end
  
  methods
    
%     function TimeStep = get.TimeStep(obj)
%       TimeStep = obj.TimeAxis(end)/obj.Length;
%     end
    
    function obj = set.TimeAxis(obj,TimeAxis)
      obj.TimeAxis = TimeAxis;
      obj = updateTimeStep(obj);
    end
  end
  
  methods(Access = private)
    
    function obj = updateTimeStep(obj)
      %Get number of points in time axis
      obj.Length = length(obj.TimeAxis);
      %Check that units are in ns
      if obj.TimeAxis(end)<50
        %If in us then rescale
        obj.TimeAxis = obj.TimeAxis/1000;
      end
      obj.TimeUnits = 'ns';
      %Get time step and round to 1ns resolution
      obj.TimeStep = round(obj.TimeAxis(end)/obj.Length,0);
    end
  end
  
end