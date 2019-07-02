classdef pdsdata
  
%==========================================================================
% Public Properties
%==========================================================================
  properties
    TimeAxis
    ExpData
    FormFactor
    DipEvoFcn
    ClusterFcn
    ModDepth
    Background
    ZeroTime
    Phase
    FitStart
  end
%==========================================================================
  
  
%==========================================================================
% Private Properties
%==========================================================================
  properties (SetAccess = private)
    TimeStep
    Length
    TimeUnits
  end
%==========================================================================
 
properties (Hidden)
    ID 
end
  
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
      ID =  java.util.UUID.randomUUID;
      obj.ID = ID.toString;
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
      if iscolumn(TimeAxis)
        TimeAxis = TimeAxis';
      end
      [obj,TimeAxis] = updateTimeStep(obj,TimeAxis);
      obj.TimeAxis = TimeAxis;
      checklengths(obj)
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
    function obj = prepare(obj,opts)
        if nargin<2 || isempty(opts)
           opts = daopts; 
        end
         obj  = prepareSignal(obj,opts);
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    function plotPreProcessing(obj)
      if ~isempty(obj.ClusterFcn)
        Figure = findobj('Tag',sprintf('ID: %s',obj.ID));
        if isempty(Figure)
          Figure = figure('Tag',sprintf('ID: %s',obj.ID),'WindowStyle','normal',...
                          'Name','Pre-processed PDS Signal');
        else
          figure(Figure);
          clf(Figure);
        end
        hold on
        plot(obj.TimeAxis,obj.ClusterFcn,'k','LineWidth',2)
        plot(obj.TimeAxis,obj.Background,'r','LineWidth',2)
        plot([0 0],[min(obj.ClusterFcn) max(obj.ClusterFcn)],'c--','LineWidth',2)
        plot([1 1]*obj.FitStart,[min(obj.ClusterFcn) max(obj.ClusterFcn)],'g--','LineWidth',2)
        axis tight
        box on, grid on
        xlabel('Time [ns]')
        ylabel('Amplitude')
        legend('Cluster Data','Background','Zero Time','FitStart')
      else
        error('Data has not been prepared.')
      end
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
    %Checks consistency of data arrays length when updating public
    %properties
    function checklengths(obj)
      if ~(isempty(obj.ExpData) || isempty(obj.TimeAxis))
        if length(obj.TimeAxis) ~=  length(obj.ExpData)
          error('ExpData and TimeAxis arrays are not equally long.')
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

