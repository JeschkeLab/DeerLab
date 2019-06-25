classdef pdsdata

%==========================================================================
  properties
    TimeAxis
    ExpData
  end
%==========================================================================

%==========================================================================
  properties (SetAccess = private)
    TimeStep
    Length
    TimeUnits
    FormFactor
    DipEvoFcn
    ClusterSignal
    ModDepth
    Background
  end
%==========================================================================  

%==========================================================================
  methods
    
    %----------------------------------------------------------------------
    function obj = set.ExpData(obj,ExpData)
      obj = getLength(obj,ExpData);
      obj.ExpData = ExpData;
      checklengths(obj)
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    function obj = set.TimeAxis(obj,TimeAxis)
      obj = getLength(obj,TimeAxis);
      obj.TimeAxis = TimeAxis;
      obj = updateTimeStep(obj);
    end
  %----------------------------------------------------------------------

  %----------------------------------------------------------------------
  function obj = prepareFormFactor(obj,Cutoff)
    %Normalize cluster signal
    obj.ClusterSignal = obj.ExpData./obj.ExpData(1);
    %Fit background
    Data2fit = obj.ClusterSignal(Cutoff:end);
    FitTimeAxis = obj.TimeAxis(Cutoff:end);
    obj.Background = fitBackground(Data2fit,obj.TimeAxis,FitTimeAxis,'exponential');
    %Correct for background by division
    obj.FormFactor = obj.ClusterSignal./obj.Background;
    obj.FormFactor = obj.FormFactor/obj.FormFactor(1);
    %Calculate modulation depth
    obj.ModDepth = 1 - obj.Background(1);
    %Get dipolar evoution function
    DipolarEvolution = obj.FormFactor - (1 - obj.ModDepth);
    DipolarEvolution = DipolarEvolution./DipolarEvolution(1);
    obj.DipEvoFcn = DipolarEvolution;
  end
  %----------------------------------------------------------------------
  
  end
%==========================================================================

%==========================================================================
  methods(Access = private)
    
    %----------------------------------------------------------------------
    function obj = updateTimeStep(obj)
      %Check that units are in ns
      if obj.TimeAxis(end)<50
        %If in us then rescale
        obj.TimeAxis = obj.TimeAxis/1000;
      end
      obj.TimeUnits = 'ns';
      %Get time step and round to 1ns resolution
      obj.TimeStep = round(obj.TimeAxis(end)/obj.Length,0);
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    function obj = getLength(obj,array)
      obj.Length = length(array);
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    %Check for equal data and time axis length
    function checklengths(obj)
      if ~(isempty(obj.ExpData) && isempty(obj.TimeAxis))
        if length(obj.TimeAxis) ~=  length(obj.ExpData)
          error('RawData and TimeAxis arrays are not equally long.')
        end
      end
    end
    %----------------------------------------------------------------------
    
  end
%==========================================================================

end