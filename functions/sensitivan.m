%
% SENSITIVAN Statistical validation of results
%
%   [median,iqr,evals] = SENSITIVAN(fcn,varpar)
%   Performs a sensibility analysis of the ouput variables returned by the
%   function (fcn) with respect to the parameter variation given in the
%   structure (varpar). The median values and inter-quartile range (median) 
%   and (iqr) of the output parameters are returned as a cell array.
%   Additionally, a third output argument (evals) can be requested, a cell 
%   array, containing the analyzed variables evaluated at each parameter 
%   combination. 
%
%   [median,iqr] = SENSITIVAN(p,vp,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order.
%
%   'AxisHandle' - Axis handle to plot the state of the validation results
%                  at each parameter permutation (default = empty)
%
%   'RandPerm' - Specifies whether to randomly permutate the validation
%   `            parameters combinations (default = true)
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function [meanOut,Upper,Lower,mainEffect,evals] = sensitivan(fcnHandle,Parameters,varargin)

if nargin<2
    error('Not enough input arguments. At least two input arguments required.')
end
if ~isa(fcnHandle,'function_handle')
    error('The first input must be a valid function handle.')
end
%Validate the required input attributes
validateattributes(Parameters,{'struct'},{'nonempty'},mfilename,'Parameters')

[AxisHandle,RandPerm] = parseoptional({'AxisHandle','RandPerm'},varargin);

validateattributes(Parameters,{'struct'},{'nonempty'},mfilename,'Parameters')


validationParam = prepvalidation(Parameters,'RandPerm',RandPerm);

% Factorial experiment design
%-----------------------------------------------------
ParNames = fieldnames(Parameters);

nout = [];
%Factorial experiment design - Start
for i=1:length(validationParam)
    %Mount input factors into user-structure
    for j=1:length(ParNames)
        argin.(ParNames{j}) =  validationParam{i,j};
    end
    %On the first run, determine the number of outputs
    if isempty(nout)
        %MATLAB's nargout returns always -1 for anonymous function handles
        nout = nargout(fcnHandle);
        if nout == -1
            notEnoughOutputs = true;
            nout = 0;
        else
            notEnoughOutputs = false;
        end
        %Iteratively increase varargout size
        while notEnoughOutputs
            try
                %If nout not large enough...
                nout = nout+1;
                varargout = cell(1,nout);
                %... this call will return an error
                [varargout{:}] = fcnHandle(argin);
                notEnoughOutputs = true;
            catch 
                nout = nout-1;
                notEnoughOutputs = false;
            end
        end
        %Prepare the output variable container
        evals = cell(1,nout);
        for ii=1:nout
            evals{i} = [];
        end
    end
    varargout = cell(1,nout);
    
    %Run the user function with current factor set
    [varargout{:}] = fcnHandle(argin);
    
    %Update statistics
    for j=1:length(varargout)
        vareval = evals{j};
        vareval(end+1,:) = varargout{j};
        %Calculate status of validation statistics
        meanOut{j} = median(vareval,1,'omitnan');
        Lower{j} = prctile(vareval,25,1);
        Upper{j} = prctile(vareval,75,1);
        evals{j} = vareval;
    end
    
    %If user passes optional plotting hook, then prepare the plot
    if ~isempty(AxisHandle)
        cla(AxisHandle)
        Ax  = 1:length(meanOut{1});
        plot(AxisHandle,Ax,meanOut{1},'k','LineWidth',1)
        hold(AxisHandle,'on')
        f = fill(AxisHandle,[Ax fliplr(Ax)] ,[Upper{1} max(fliplr(Lower{1}),0)],...
            'b','LineStyle','none');
        f.FaceAlpha = 0.5;
        hold(AxisHandle,'off')
        axis(AxisHandle,'tight')
        grid(AxisHandle,'on')
        box(AxisHandle,'on')
        title(sprintf('Run %i/%i',i,length(validationParam)))
        drawnow
    end
end
%Factorial experiment design - End


%Main effect analysis
%-----------------------------------------------------
for i = 1:nout
    data = evals{i};
    for j = 1:size(validationParam,2)
        clear evalmean
        clear set
        subset = validationParam(1:size(data,1),j);
        if ischar(subset{1})
            uni = unique(subset);
        else
            subset = cell2mat(subset);
            uni = unique(subset);
        end
        for ii=1:length(uni)
            if iscell(subset)
                idx =  find(contains(subset,uni{ii}));
            else
                idx = find(subset==uni(ii));
            end
            evalmean(ii,:) = pdist(data(idx,:),'euclidean');
        end
        if size(evalmean,1)<=2
            mainEffect{i}(j) = abs(mean((evalmean(1,:) - evalmean(2,:))));
        else
            mainEffect{i}(j) = abs(mean(((evalmean(1,:) - evalmean(2,:)) + (evalmean(1,:)) - evalmean(3,:)) + (evalmean(2,:) - evalmean(3,:))));
        end
    end
end


if nout==1
    meanOut = meanOut{1};
    Upper = Upper{1};
    Lower = Lower{1};
    evals = evals{1};
    mainEffect = mainEffect{1};
end








