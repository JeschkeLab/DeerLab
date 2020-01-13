%
% SENSITIVAN Sensitivity analysis by factorial design
%
%   [med,up,lo] = SENSITIVAN(fcn,varpar)
%   Performs a sensibility analysis of the ouput variables returned by the
%   function (fcn) with respect to the parameter variation given in the
%   structure (varpar). The median values (med), as well as the lower (lo)
%   and upper (up) values of the output parameters are returned as a cell array.
%
%   [~,~,~,main,inter] = SENSITIVAN(fcn,varpar)
%   Additionally, the main effects (main) and interaction (inter) between
%   the factors in varpar are returned.
%
%   [~,~,~,~,~,evals] = SENSITIVAN(fcn,varpar)
%   The last output argument returns a cell array containing all of the
%   outputs computed for all factor level combinations.
%
%   [args] = SENSITIVAN(fcn,varpar,'Property',Value)
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


function [meanOut,Upper,Lower,mainEffect,Interaction,evals] = sensitivan(fcnHandle,Parameters,varargin)

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
%Get names of the variables givne by the user
ParNames = fieldnames(Parameters);

nout = [];
%Factorial experiment design - Start
for i=1:size(validationParam,1)
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
        %Iteratively increase varargout size to determine number of outputs
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
        meanOut = cell(1,nout);
        Upper = cell(1,nout);
        Lower = cell(1,nout);
        for ii=1:nout
            evals{i} = [];
        end
    end
    varargout = cell(1,nout);
    
    %Run the user function with current factor set
    [varargout{:}] = fcnHandle(argin);
    
    %Check that the outputs are strictly numerical
    anycharoutput = cellfun(@(x)ischar(x),varargout);
    if anycharoutput
        error('String output arguments f by the input function handle are not accepted by sensitivan.')
    end
    
    %Since MATLAB does not know the size of the varargout variables...
    for j=1:length(varargout)
        %... extract them one by one
        vareval = evals{j};
        out = varargout{j};
        %Format all vectors as columns
        if isvector(out) && iscolumn(out)
            out = out.';
        end
        %Check the consistency of the output variabe size
        if ~isempty(vareval)
            sizerror = false;
           if isvector(out) && numel(out)~=numel(vareval(1,:))
               sizeerror = true;
           elseif ismatrix(out) && any(size(out)~=size(squeeze(vareval(1,:,:,:,:))))
               sizeerror = true;
           end
           %If size is inconsistent, launch error with advice
           if sizerror
              error(['Inconsistent output variable size. ',...
                    'One of the outputs of the analyzed function is changing its size in between runs. ',...
                    'To solve this, fix the axis of the output and interpolate the result. \n%s'],...
                    '  Ex/ outFix = interp1(varAxis,out,fixAxis,''spline'')') 
           end
        end
        %... and store them in a N-dimensional container
        %The unused singlet dimensions are automatically ignored by MATLAB
        vareval(end+1,:,:,:,:) = out;
        evals{j} = vareval;
    end
    
    %Once at least two runs have been completed (required for statistics)
    if i>1
        %Update statistics
        for j=1:length(varargout)
            vareval = evals{j};
            %Calculate status of sensitivity analysis statistics
            meanOut{j} = squeeze(median(vareval,1,'omitnan'));
            Lower{j} = percentile(vareval,25,1).';
            Upper{j} = percentile(vareval,75,1).';
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
    %Finished, proceed to next combination of factors
end
%Factorial experiment design - End


% Factors main effect analysis
%-----------------------------------------------------
%Run through all function output variables
for i = 1:nout
    %Get all evaluations of that variable
    data = evals{i};
    %Loop through all factors
    for j = 1:size(validationParam,2)
        clear evalmean
        clear set
        subset = validationParam(1:size(data,1),j);
        %Depending on data type, find unique factor levels
        if isa(subset{1},'function_handle')
            subset = cellfun(@func2str,subset,'UniformOutput',false);
        end
        if ischar(subset{1})
            uni = unique(subset);
        else
            subset = cell2mat(subset);
            uni = unique(subset);
        end
        %Loop throught the levels of the factor
        for ii=1:length(uni)
            
            %Identify the indices of the evaluations using that level
            if iscell(subset)
                idx =  find(contains(subset,uni{ii}));
            else
                idx = find(subset==uni(ii));
            end
            %Get the subset of data for that factor level
            M = data(idx,:).';
            %Construct a Euclidean distance map (kind of a toolbox-free pdist)
            map = bsxfun(@plus,dot(M,M,1),dot(M,M,1)')-2*(M'*M);
            %Make it an upper triangle-matrix...
            map = triu(map,1);
            %...and get the mean Euclidean distance
            evalmean(ii) = mean(map(map~=0));
        end
        %The main effect is given by the difference between mean Euclidean
        %distances at different factor levels
        main{i}{j} = abs(diff(evalmean));
    end
end

%Format the output into a nice structure with the names given by the user
for i=1:length(main)
    for j=1:length(ParNames)
        mainEffect{i}.(ParNames{j}) =  main{i}{j};
    end
end

% Factors interaction analysis (beta)
%-----------------------------------------------------
for i = 1:nout
    data = evals{i};
    for j = 1:size(validationParam,2)
        for k = 1:size(validationParam,2)
            clear evalmean
            clear set
            clear subset
            clear uni
            clear main
            %Get subsets of the two interacting factors
            subset{1} = validationParam(1:size(data,1),j);
            subset{2} = validationParam(1:size(data,1),k);
            %Get the primary factor main effects for both levels of the
            %secondary factor
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
                %Get main effects for upper and lower levels
                main(jj) = abs(evalmean(1) - evalmean(2));
                
            end
            
            %Compute the interaction between the factors
            Interaction{i}(j,k) = abs(main(1) - main(2));
            
        end
    end
end

%If only one output variable has been evaluated, then don't return cell array
if nout==1
    meanOut = meanOut{1};
    Upper = Upper{1};
    Lower = Lower{1};
    evals = evals{1};
    mainEffect = mainEffect{1};
    Interaction = Interaction{1};
end


end

%Construct free workaround for prctile function (Statistics Toolbox)
function pct = percentile(M,p,dim)

%Set requested dimension as the first dimension
dimIdx = 1:ndims(M);
dimIdx = dimIdx(dimIdx~=dim);
M = permute(M,[dim dimIdx]);

sizeM = size(M);

%Vectorize all other dimensions
if numel(sizeM)>2
    V = reshape(M,[sizeM(1),prod(sizeM(2:end))]);
else
    V = M;
end

%Prepare prctile function to compute throughout first dimension
prctile = @(v,p) interp1(linspace(0.5/length(v), 1-0.5/length(v), length(v))', sort(v), p*0.01, 'pchip');
colfcn = @(func, matrix) @(col) func(matrix(:,col));
colfcn = @(func, matrix) arrayfun(colfcn(func, matrix), 1:size(matrix,2), 'UniformOutput', false)';
takeall = @(x) reshape([x{:}], size(x{1},2), size(x,1))';
applyrowfcn = @(func, matrix) takeall(colfcn(func, matrix));

%Apply function nest to vectorized matrix
pct = applyrowfcn(@(M)prctile(M,p),V);

if numel(sizeM)>2
    %Reshape results back to original size
    pct = reshape(pct,sizeM(2:end));
end

end





