classdef distdistr
    
   properties 
      Distribution
      DistanceAxis
      ExpSignal
      FitSignal
      TimeAxis
   end
   
   properties (Access = private)
      Method
   end
   
   properties (Hidden)
      ID 
   end
    
   methods
        function obj = distdistr(varargin)
            if nargin>0
                if mod(nargin,2)
                    error('Wrong number of arguments')
                end
                for i=1:2:length(varargin)
                    obj = setProperty(obj,varargin{i},varargin{i+1});
                end
            end
            JavaID =  java.util.UUID.randomUUID;
            obj.ID = JavaID.toString;
        end 
        
        function obj = set.Distribution(obj,Distribution)
            CallStack = dbstack(1);
            obj.Method = CallStack(1).name;
            obj.Distribution = Distribution;
        end
        
    function plotDistr(obj)
      if ~isempty(obj.Distribution)
        Figure = findobj('Tag',sprintf('ID: %s',obj.ID));
        if isempty(Figure)
          Figure = figure('Tag',sprintf('ID: %s',obj.ID),'WindowStyle','normal',...
                          'Name','Distance Distribution');
        else
          figure(Figure);
          clf(Figure);
        end
        subplot(121)
        plot(obj.DistanceAxis,obj.Distribution,'k','LineWidth',2)
        axis tight
        xlabel('Distance [nm]')
        ylabel('P(r)')
        subplot(122)
        hold on
        plot(obj.TimeAxis,obj.ExpSignal,'k','LineWidth',2)
        plot(obj.TimeAxis,obj.FitSignal,'r','LineWidth',2)
        axis tight
        box on, grid on
        xlabel('Time [ns]')
        ylabel('Amplitude')
        legend('Input','Fit')
      else
        error('Data has not been prepared.')
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