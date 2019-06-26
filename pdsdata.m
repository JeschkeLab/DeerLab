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
    ClusterFcn
    ModDepth
    Background
  end
  %==========================================================================
  
  %==========================================================================
  methods
    
    %----------------------------------------------------------------------
    function obj = pdsdata(varargin)
      if nargin>0
        if mod(nargin,2)
          error('Wrong number of arguments')
        end
        for i=1:2:length(varargin)
          obj = setProperty(obj,varargin{i},varargin{i+1});
        end
      end
    end
    %----------------------------------------------------------------------
    
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
      Cutoff = 70;
      %Normalize cluster signal
      obj.ClusterFcn = obj.ExpData/obj.ExpData(1);
      
      %Fit background
      Data2fit = obj.ClusterFcn(Cutoff:end);
      FitTimeAxis = obj.TimeAxis(Cutoff:end);
      obj.Background = fitBackground(Data2fit,obj.TimeAxis,FitTimeAxis,'exponential');
      
      %Correct for background by division
      obj.FormFactor = obj.ClusterFcn./obj.Background;
      obj.FormFactor = obj.FormFactor/obj.FormFactor(1);
      
      %Calculate modulation depth
      obj.ModDepth = 1 - obj.Background(1);
      
      %Get dipolar evoution function
      DipolarEvolution = obj.FormFactor - (1 - obj.ModDepth);
      DipolarEvolution = DipolarEvolution./DipolarEvolution(1);
      obj.DipEvoFcn = DipolarEvolution;
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    function obj = plot(obj)
      obj = getLength(obj,TimeAxis);
      obj.TimeAxis = TimeAxis;
      obj = updateTimeStep(obj);
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
    
    %----------------------------------------------------------------------
    function obj = setProperty(obj,Property,Value)
      if ~isa(Property,'char') && ~isa(Property,'string')
        error('Invalid first data argument')
      end
      eval(sprintf('obj.%s = Value;',Property));
    end
    %----------------------------------------------------------------------
    
    
  end
  %==========================================================================
  
end