%
% VALIDATE Statistical validation of results
%
%   [mean,std,evals] = VALIDATE(fcn,varpar)
%   Performs a sensibility analysis of the ouput variables returned by the
%   function (fcn) with respect to the parameter variation given in the
%   structure (varpar). The mean values and standard deviations (mean) and 
%   (std) of the output parameters are returned as a cell array.
%   Additionally, a output third argument (evals) can be requested, a cell 
%   array, containing the analyzed variables evaluated at each parameter 
%   combination. 
%
%   [mean,std] = VALIDATE(p,vp,'Property',Value)
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
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.


function [meanOut,stdOut,evals] = validate(fcnHandle,Parameters,varargin)

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

%-----------------------------------------------------
% Run validation
%-----------------------------------------------------

ParNames = fieldnames(Parameters);

%Prepare the output variable container
evals = cell(1,length(validationParam));
for i=1:length(validationParam)
    evals{i} = [];
end
nout = [];
%Run over all validation parameter permutations
for i=1:length(validationParam)
    for j=1:length(ParNames)
        argin.(ParNames{j}) =  validationParam{i,j};
    end
    if isempty(nout)
        nout = nargout(fcnHandle);
        if nout == -1
            notEnoughOutputs = true;
            nout = 1;
        else
            notEnoughOutputs = false;
        end
        while notEnoughOutputs
            try
                nout = nout+1;
                varargout = cell(1,nout);
                %Run the build function
                [varargout{:}] = fcnHandle(argin);
                notEnoughOutputs = false;
            catch
                notEnoughOutputs = true;
            end
        end
    end
    varargout = cell(1,nout);
    [varargout{:}] = fcnHandle(argin);
    for j=1:length(varargout)
        vareval = evals{j};
        vareval(end+1,:) = varargout{j};
        %Calculate status of validation statistics
        meanOut{j} = mean(vareval,1);
        stdOut{j} = std(vareval,[],1);
        evals{j} = vareval;
    end
    %If user passes optional plotting hook, then prepare the plot
    if ~isempty(AxisHandle)
        cla(AxisHandle)
        Ax  = 1:length(meanOut{1});
        plot(AxisHandle,Ax,meanOut{1},'k','LineWidth',1)
        hold(AxisHandle,'on')
        f = fill(AxisHandle,[Ax fliplr(Ax)] ,[meanOut{1}+stdOut{1} max(fliplr(meanOut{1}-stdOut{1}),0)],...
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

if nout==1
    meanOut = meanOut{1}.';
    stdOut = stdOut{1}.';
    evals = evals{1};
end








