classdef pdsdata
  
%==========================================================================
% Public Properties
%==========================================================================
  properties
    TimeAxis
    ExpData
  end
%==========================================================================
  
  
%==========================================================================
% Private Properties
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
    ZeroTime
    Phase
  end
%==========================================================================
  
  
%==========================================================================
% Public methods
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
      if iscolumn(ExpData)
        ExpData = ExpData';
      end
      obj.ExpData = ExpData;
      checklengths(obj)
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    function obj = set.TimeAxis(obj,TimeAxis)
      obj = getLength(obj,TimeAxis);
      [obj,TimeAxis] = updateTimeStep(obj,TimeAxis);
      obj.TimeAxis = TimeAxis;
    end
    %----------------------------------------------------------------------
  
    %----------------------------------------------------------------------
    function obj = set.ModDepth(obj,ModDepth)
      if ModDepth<0 || ModDepth>1
        error('Modulation depth is out of boundaries.')
      end
      obj.ModDepth = ModDepth;
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    function obj = set.Length(obj,Length)
      if Length==0
        error('Input array is empty.')
      end
      obj.Length = Length;
    end
    %----------------------------------------------------------------------

    %----------------------------------------------------------------------
    function obj = prepareFormFactor(obj)
      if isempty(obj.ExpData)
        error('ExpData property is empty.')
      end
      if isempty(obj.TimeAxis)
        error('TimeAxis property is empty.')
      end
      Cutoff = 70;
      %Normalize cluster signal
      obj.ClusterFcn = obj.ExpData/obj.ExpData(1);
      [obj.ClusterFcn,obj.Phase] = correctPhase(obj.ClusterFcn);
      [obj.ClusterFcn,obj.TimeAxis,obj.ZeroTime] = correctZeroTime(obj.ClusterFcn,obj.TimeAxis);
       
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
    function plot(obj)
        plot(obj.TimeAxis,obj.ClusterFcn)
    end
    %----------------------------------------------------------------------
  
  end
%==========================================================================
  
%==========================================================================
% Private methods
%==========================================================================
  methods(Access = private)
    
    %----------------------------------------------------------------------
    function [obj,TimeAxis] = updateTimeStep(obj,TimeAxis)
      %Check that units are in ns
      if TimeAxis(end)<50
        %If in us then rescale
        TimeAxis = TimeAxis*1000;
      end
      obj.TimeUnits = 'ns';
      %Get time step and round to 1ns resolution
      obj.TimeStep = round(TimeAxis(end)/obj.Length,0);
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

