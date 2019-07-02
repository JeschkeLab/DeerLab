classdef distdistr
    
   properties 
      Distribution
      DistanceAxis
      ExpSignal
      TimeAxis
      RegParam
   end
   
   properties (SetAccess = private)
      Method
      FitSignal
      Length
   end
   
   properties (Hidden)
      ID 
      SignalID
   end
    
   methods
        function obj = distdistr(varargin)
            %Check if user given inputs
            if nargin>0
                %Only pairwise inputs accepted
                if mod(nargin,2)
                    error('Wrong number of arguments')
                end
                %Set properties given to constructor
                for i=1:2:length(varargin)
                    obj = setProperty(obj,varargin{i},varargin{i+1});
                end
            end
            %Give unique ID to the variable
            JavaID =  java.util.UUID.randomUUID;
            obj.ID = JavaID.toString;
        end 
        
        function obj = set.Distribution(obj,Distribution)
            %Get function call stack
            CallStack = dbstack;
            %Look for function calling the class functions
            i = 1;
            while ~isempty(strfind(CallStack(i).name,'distdistr'))
            i  = i + 1;
            end
            %Store which function/method constructed the distribution
            obj.Method = CallStack(i).name;
            %Check distribution array
            if ~iscolumn(Distribution)
               Distribution = Distribution'; 
            end
            obj.Length = length(Distribution);
            %Force integral normalization
            Distribution = Distribution/sum(Distribution);
            obj.Distribution = Distribution;
        end
                
    function plot(obj,CompDistAxis,CompDistribution)
      if ~isempty(obj.Distribution)
        Figure = findobj('Tag',sprintf('ID: %s',obj.ID));
        if isempty(Figure)
          Figure = figure('Tag',sprintf('ID: %s',obj.ID),'WindowStyle','normal',...
                          'Name','Distance Distribution');
        else
          figure(Figure);
          clf(Figure);
        end
        %Plot distribution
        subplot(121)
        plot(obj.DistanceAxis,obj.Distribution,'r','LineWidth',2)
        %If requested plot against input distribution
        if nargin>1
        hold on
        plot(CompDistAxis,CompDistribution,'k','LineWidth',2)
        end
        axis tight, box on, grid on
        xlabel('Distance [nm]'),ylabel('P(r)')
        %Plot 
        subplot(122)
        hold on
        plot(obj.TimeAxis,obj.ExpSignal,'k','LineWidth',2)
        plot(obj.TimeAxis,obj.FitSignal,'r','LineWidth',2)
        axis tight, box on, grid on
        xlabel('Time [ns]'),ylabel('Amplitude')
        legend('Input','Fit')
      else
        error('Data has not been prepared.')
      end
    end

    function obj = getFitSignal(obj)
        TimeStep = obj.TimeAxis(2) - obj.TimeAxis(1);
        Kernel = getKernel(obj.Length,TimeStep);
        obj.FitSignal = Kernel*obj.Distribution;
    end
    
    function match = findSignal(obj) 
        baseVariables = evalin('base' , 'whos');
        match = [];
        for i = 1:length(baseVariables)
            if (strcmpi(baseVariables(i).class , 'pdsdata')) % compare classnames
                variable = evalin('base',baseVariables(i).name);
                querySignalID = getfield(variable,'ID');
                if obj.SignalID == querySignalID
                    match = variable;
                end
            end
        end
        if isempty(match)
           error('Corresponding pdsdata class object not found or deleted.') 
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