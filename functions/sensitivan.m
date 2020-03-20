%
% SENSITIVAN Sensitivity analysis by factorial design
%
%   [stats] = SENSITIVAN(fcn,varpar)
%   [stats,factors] = SENSITIVAN(...)
%   [stats,factors,evals] = SENSITIVAN(...)
%   ... = SENSITIVAN(fcn,varpar,'Property',Value)
%
%   Performs a sensitivity and factor analysis of the output variables returned
%   by the function (fcn) with respect to the parameter variation given in the
%   structure (varpar).
%
%   Inputs:
%     fcn     function handle with the function to call
%     varpar  structure with parameters to vary
%
%     Property-value pairs (optional)
%     'AxisHandle' - Axis handle to plot the state of the validation results
%                    at each parameter permutation. If empty, nothing is plotted.
%                   (default = empty)
%
%     'RandPerm' - Specifies whether to randomly permutate the validation
%                  parameters combinations (default = true)
%   
%     'dynamicStats' - Specifies whether the statistical estimators are
%                      computed using the full set of observations or approximated
%                      dynamically at each iteration. 
%                      (default = depends on available memory)
%
%   Outputs:
%     stats     structure array with summary statistics for each variable
%       .median medians of the output variables
%       .mean   means of the output variables
%       .std    standard deviations of the output variables
%       .p2     2nd  percentiles of the output variables
%       .p25    25th percentiles of the output variables
%       .p75    75th percentiles of the output variables
%       .p98    98th percentiles of the output variables

%     factors   results of factor analysis
%       .main   main effects
%       .inter  interactions between factors
%
%     evals     cell array with all outputs computed for all factor level
%               combinations
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function [stats,factors,evals] = sensitivan(fcnHandle,Parameters,varargin)

% Input validation
%-------------------------------------------------------------------------------
if nargin<2
    error('Not enough input arguments. At least two input arguments required.');
end
if ~isa(fcnHandle,'function_handle')
    error('The first input must be a function handle.');
end
try
    nargin(fcnHandle);
catch
    error('The function given in the first input cannot be found.');
end

validateattributes(Parameters,{'struct'},{'nonempty'},mfilename,'Parameters');

[AxisHandle,RandPerm,dynamicStats] = ...
    parseoptional({'AxisHandle','RandPerm','dynamicStats'},varargin);

if ~isempty(dynamicStats)
    validateattributes(dynamicStats,{'logical'},{'nonempty'},mfilename,'dynamicStats')
end

% Get parameter names, generate all parameter combinations
ParNames = fieldnames(Parameters);
nParameters = numel(ParNames);

ParamList = prepvalidation(Parameters,'RandPerm',RandPerm);
nCombinations = size(ParamList,1);

% Factorial experiment design
%-------------------------------------------------------------------------------
% Get names of the variables given by the user

counter = 0;
nout = [];

for i = 1:nCombinations
    
    % Assemble input factors into user structure
    for p = 1:nParameters
        argin.(ParNames{p}) =  ParamList{i,p};
    end
    
    % On the first run, determine the number of outputs and allocate arrays
    if isempty(nout)
        nout = getmaxnargout(fcnHandle,argin);
        %Pre-allocate memory for stats structure array
        stats = repmat(struct('median',[],'mean',[],'std',[],'p2',[],'p25',[],'p75',[],'p98',[]),nout,1);
        evals = cell(1,nout);
        sizeOut = cell(1,nout);
        %Prepare containers for dynamic statistical estimators
        markers2 = cell(1,nout);
        markers25 = cell(1,nout);
        markers50 = cell(1,nout);
        markers75 = cell(1,nout);
        markers98 = cell(1,nout);
        OutSum(1:nout) = {0};
        OutSumSq(1:nout) = {0};
        
    end
    
    % Run the user function with current factor set
    varargsout = cell(1,nout);
    [varargsout{:}] = fcnHandle(argin);
    
    % Assert that all outputs are strictly numerical
    numericOutput = cellfun(@(x)isnumeric(x),varargsout);
    if ~all(numericOutput)
        error('Non-numeric output arguments by the input function are not accepted.');
    end
    
    %Perform memory checks to estimate if there is enough memory
    if i==1
        isEnoughMemory = memorycheck(varargsout,nargout,nCombinations);
    end
    
    %If not specified by the user, decide which statistics estimators to use
    if isempty(dynamicStats) && isEnoughMemory
        dynamicStats = false;
    elseif isempty(dynamicStats) && ~isEnoughMemory
        dynamicStats = true;
    elseif dynamicStats && ~isEnoughMemory
        error(['The current sensitivity analysis using order statistics ' ...
               'requires more memory than available. Change the option ' ...
               '''DynamicStats'' to true to proceed.'])
    end
    
    % Store outputs in an N-dimensional array
    for iOut = 1:nout
        out = varargsout{iOut};
        
        if isempty(sizeOut{iOut})
            sizeOut{iOut} = size(out);
        elseif any(sizeOut{iOut}~= size(out))
            error(['Inconsistent output variable size. ',...
                'One of the outputs of the analyzed function is changing its size in between runs. ',...
                'To solve this, fix the axis of the output and interpolate the result. \n%s'],...
                '  Ex/ outFix = interp1(varAxis,out,fixAxis,''pchip'')')
        end
        % Convert vectors to rows
        if iscolumn(out)
            out = out.';
        end
        

        if isEnoughMemory && ~dynamicStats
            %Store the outputs of the current parameter combination
            evals{iOut} = cat(1,evals{iOut},shiftdim(out,-1));
        end
    end
    
    % Update statistics
    %---------------------------------------------------------------
    counter = counter + 1;
    % Evalutate the costly percentile function only every 5 combinations or
    % after all combinations have been evaluated
    if counter == 5 || i == nCombinations
        
        if dynamicStats
            %Dynamically update the statistical estimators for all outputs
            for j = 1:nout
                %Get current sample/observation for output j
                sample = varargsout{j};
                
                %Dynamic mean and standard deviation
                OutSum{j} = OutSum{j} + sample;
                OutSumSq{j} = OutSumSq{j} + sample.*sample;
                stats(j).mean = OutSum{j}/i;
                stats(j).std = OutSumSq{j} - (OutSum{j}.*OutSum{j}/i)/(i-1);
                
                %Dynamic percentiles/quantiles
                [stats(j).median, markers50{j}] = dynprctile(0.50,sample,markers50{j});
                [stats(j).p98, markers98{j}] = dynprctile(0.98,sample,markers98{j});
                [stats(j).p75, markers75{j}] = dynprctile(0.75,sample,markers75{j});
                [stats(j).p25, markers25{j}] = dynprctile(0.25,sample,markers25{j});
                [stats(j).p2, markers2{j}] = dynprctile(0.02,sample,markers2{j});
            end
        else
            %Update statistical estimators using order statistics
            for j = 1:nout
                vareval = evals{j};
                stats(j).median = squeeze(median(vareval,1,'omitnan'));
                stats(j).mean = squeeze(mean(vareval,1,'omitnan'));
                stats(j).std = squeeze(std(vareval,0,1,'omitnan'));
                if i>1
                    stats(j).p2  = percentile(vareval,2,1).';
                    stats(j).p25 = percentile(vareval,25,1).';
                    stats(j).p75 = percentile(vareval,75,1).';
                    stats(j).p98 = percentile(vareval,98,1).';
                end
            end
        end
        %Reset counter
        counter = 0;
    end
    
    % If user passes optional plotting hook, then prepare the plot
    if ~isempty(AxisHandle) && i>1
        cla(AxisHandle)
        Ax = 1:length(stats(1).median);
        plot(AxisHandle,Ax,stats(1).median,'k','LineWidth',1)
        hold(AxisHandle,'on')
        f = fill(AxisHandle,[Ax fliplr(Ax)] ,[stats(1).p75; max(flipud(stats(1).p25),0)],...
            'b','LineStyle','none');
        f.FaceAlpha = 0.5;
        hold(AxisHandle,'off')
        axis(AxisHandle,'tight')
        grid(AxisHandle,'on')
        box(AxisHandle,'on')
        title(sprintf('Run %i/%i',i,length(ParamList)));
        drawnow
    end
    
end


%If requested, proceed with factor analysis
if nargout>1
    % Factors main effect analysis
    %-------------------------------------------------------------------------------
    main(1:nout) = {cell(nParameters,1)};
    mainEffect(1:nout) = {struct()};
    % Loop over all function output variables
    for iOut = 1:nout
        % Get all evaluations of that variable
        data = evals{iOut};
        % Loop over all factors
        for iPar = 1:nParameters
            clear evalmean set
            subset = ParamList(1:size(data,1),iPar);
            % Find unique factor levels
            if isa(subset{1},'function_handle')
                subset = cellfun(@func2str,subset,'UniformOutput',false);
            elseif ~ischar(subset{1})
                subset = cell2mat(subset);
            end
            uni = unique(subset);
            
            % Loop throught the levels of the factor
            evalmean = zeros(numel(uni),1);
            for ii = 1:numel(uni)
                
                % Identify the indices of the evaluations using that level
                if iscell(subset)
                    idx =  find(contains(subset,uni{ii}));
                else
                    idx = find(subset==uni(ii));
                end
                % Get the subset of data for that factor level
                M = data(idx,:).';
                % Construct a Euclidean distance map
                map = bsxfun(@plus,dot(M,M,1),dot(M,M,1)')-2*(M'*M);
                % Make it an upper triangle-matrix and get the mean Euclidean distance
                map = triu(map,1);
                evalmean(ii) = mean(map(map~=0));
            end
            % The main effect is given by the difference between mean Euclidean
            % distances at different factor levels
            main{iOut}{iPar} = abs(diff(evalmean));
        end
    end
    
    % Convert cell array to strucure array with original parameter names
    for i = 1:nout
        for p = 1:nParameters
            mainEffect{i}.(ParNames{p}) =  main{i}{p};
        end
    end
    factors.main = mainEffect;
    
    
    % Factors interaction analysis
    %-------------------------------------------------------------------------------
    Interaction(1:nout) = {zeros(nParameters)};
    for i = 1:nout
        data = evals{i};
        for j = 1:nParameters
            for k = 1:nParameters
                clear evalmean set subset uni main
                % Get subsets of the two interacting factors
                subset{1} = ParamList(1:size(data,1),j);
                subset{2} = ParamList(1:size(data,1),k);
                % Get the primary factor main effects for both levels of the
                % secondary factor
                for jj=1:length(subset)
                    tmp = subset{jj};
                    if isa(tmp{1},'function_handle')
                        tmp = cellfun(@func2str,tmp,'UniformOutput',false);
                    end
                    if ischar(tmp{1})
                        uni{jj} = unique(tmp);
                    else
                        tmp = cell2mat(tmp);
                        uni{jj} = unique(tmp);
                    end
                    subset{jj} = tmp;
                    
                    for ii=1:length(uni{jj})
                        unitmp = uni{jj};
                        if iscell(subset{jj})
                            idx =  find(contains(subset{jj},unitmp{ii}));
                        else
                            idx = find(subset{jj}==unitmp(ii));
                        end
                        M = data(idx,:).';
                        map = bsxfun(@plus,dot(M,M,1),dot(M,M,1)')-2*(M'*M);
                        map = triu(map,1);
                        evalmean(ii) = mean(map(map~=0));
                    end
                    % Get main effects for upper and lower levels
                    main{jj} = abs(diff(evalmean));
                end
                
                % Compute the interaction between the factors
                Interaction{i}(j,k) = abs(mean(main{1}) - mean(main{2}));
                
            end
        end
    end
    factors.inter = Interaction;
else
    factors.inter = cell(1,1);
    factors.main = cell(1,1);
end

% If only one output variable has been evaluated, then don't return cell array
if nout==1
    factors.main = factors.main{1};
    factors.inter = factors.inter{1};
    evals = evals{1};
end

end

function nout = getmaxnargout(fcnHandle,argin)
nout = nargout(fcnHandle);
variableOutputs = nout<0;

if ~variableOutputs
    return;
end

% If the functions defines a variable number of outputs, iteratively increase
% number of requested outputs until the function crashes, to determine maximum
% number of outputs.
nout = abs(nout)-1;
done = false;
while ~done
    try
        nout = nout+1;
        varargout = cell(1,nout);
        [varargout{:}] = fcnHandle(argin);
    catch
        nout = nout-1;
        done = true;
    end
end
end

% Calculate percentile (similar to prctile function in Statistics Toolbox)
function Y = percentile(X,p,dim)

% Set requested dimension as the first dimension
dimIdx = 1:ndims(X);
dimIdx = dimIdx(dimIdx~=dim);
X = permute(X,[dim dimIdx]);

% Get size of data
sizeX = size(X);

% Vectorize all other dimensions
if numel(sizeX)>2
    X = reshape(X,[sizeX(1),prod(sizeX(2:end))]);
end

N = size(X,1);
% Sort data about first dimension
X = sort(X,1);
% Get list of available percentiles
pList = 100*(0.5:1:N-0.5)/N;
% Interpolate from list to requested percentile
Y = interp1(pList,X,p,'linear','extrap');

if numel(sizeX)>2
    % Reshape results back to original size
    Y = reshape(Y,sizeX(2:end));
end

end

function isEnoughMemory = memorycheck(sample,nargout,nCombinations)

%Estimate the memory costs of the current sensitivity analysis
info = whos('sample');
memoryusage = nCombinations*info.bytes/1e9; %GB
[memoryfree,~] = memorystatus(); %GB

%Stop now if the estimate would exceed virtual memory available
%to MATLAB. Ootherwise it will crash mid-execution.
if memoryusage > memoryfree && nargout>1
    error(['The current sensitivity analysis is requesting a factor '...
        'analysis and/or the evaluated outputs at all combinations.'...
        ' According to the input this will require %.2f GB of memory, '...
        'which exceed the amoun available to MATLAB (%.2f GB).'],round(memoryusage/1e9,2),round(memoryfree/1e9,2))
elseif  memoryusage > memoryfree && nargout==1
    isEnoughMemory = false;
else
    isEnoughMemory = true;
end

%If the estimate exceeds 3GB, use dynamic
if memoryusage/1e9 > 3 && nargout>1
    warning('on','all')
    warning(['sensitivan will require %.2f GB of memory. This may slow down '...
        'your analysis or MATLAB instance. Consider reducing the amount '...
        'of factors/levels or requesting only the ''stats'' output.'],memoryusage/1e9)
    warning('off','all')
elseif memoryusage/1e9 > 3 && nargout==1
    isEnoughMemory = false;
end
end
