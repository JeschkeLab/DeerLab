classdef distribution < DeerAnalysis_object
    
   properties 
      Distribution
      DistanceAxis
   end
   
   properties (SetAccess = private)
      FitSignal
      Length
   end
   
   properties (Hidden)
      ID 
   end
    
   methods
        function obj = distribution(varargin)           
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
        
        function value = get.Length(obj)
            value = length(obj.Distribution);
        end
        
        function obj = set.Distribution(obj,Distribution)
            %Check distribution array
            if ~iscolumn(Distribution)
               Distribution = Distribution'; 
            end
            Distribution = normalize(Distribution);
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
        legendTag = {};
        %If requested plot against input distribution
        if nargin>1
            hold on
            plot(CompDistAxis,CompDistribution,'k','LineWidth',2)
            legendTag{1} = 'Model';
        end
        plot(obj.DistanceAxis,obj.Distribution,'r','LineWidth',2)
        axis tight, box on, grid on
        xlabel('Distance [nm]'),ylabel('P(r) [nm^{-1}]')
        legendTag{end+1} = 'Fit';
        legend(legendTag)
        %Plot 
        subplot(122)
        hold on
        plot(obj.TimeAxis,obj.Signal,'k','LineWidth',2)
        plot(obj.TimeAxis,obj.FitSignal,'r','LineWidth',2)
        axis tight, box on, grid on
        xlabel('Time [ns]'),ylabel('Amplitude')
        legend('Input','Fit')
      else
        error('Data has not been prepared.')
      end
    end

    function Distribution = normalize(Distribution)
        Distribution = Distribution/sum(Distribution);
    end
    
%     function FitSignal = get.FitSignal(obj)
%         Kernel = queryKernel(obj);
%         FitSignal = Kernel*obj.Distribution;
%     end
%         function Kernel = queryKernel(obj)
%         Kernel = queryKernel@dataset;
%     end
    
   end
   
   methods (Static)


    
   end
 
    
end