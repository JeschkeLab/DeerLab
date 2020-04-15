function stats = bootan(fcn,Vexp,Vfit,varargin)

% Parse the optional parameters in varargin
optionalProperties = {'Verbose','Samples'};
[Verbose,nSamples] = parseoptional(optionalProperties,varargin);

if isempty(nSamples)
    nSamples = 1000;
else
    validateattributes(nSamples,{'numeric'},{'nonnegative','scalar','nonempty'},mfilename,'Nsamples')
end

if isempty(Verbose)
    Verbose = false;
else
    validateattributes(Verbose,{'logical'},{'nonempty'},mfilename,'Verbose')
end

if numel(Vexp)~=numel(Vfit)
    error('V and Vfit (2nd and 3rd input) must have the same number of element.')
end

% Determine the number of outputs returned by function handle and allocate arrays
nout = getmaxnargout(fcn,Vexp);
evals = cell(1,nout);
sizeOut = cell(1,nout);

% Get residuals and estimate standard deviation
residuals = Vfit - Vexp;
sigma = std(residuals);

% Generate bootstrap samples from data
%-------------------------------------------------------------------------------
warning('off','all')
for iSample = 1:nSamples
    
    % Inform of progress if requested
    if Verbose
        S = sprintf('Bootstrapping: %i/%i samples finished',iSample,nSamples);
        fprintf(S);
    end
    
    % Get a bootstrap sample
    if iSample>1
        % Resample residuals (using the fact that they are Gaussian distributed
        % and use data residuals to get variance)
        Vresample = Vfit + sigma*randn(size(Vfit));
    else
        % Use original data
        Vresample = Vexp;
    end
    
    % Run the model function with bootstrap sample
    varargsout = cell(1,nout);
    [varargsout{:}] = fcn(Vresample);
    
    % Assert that all outputs are strictly numerical
    numericOutput = cellfun(@(x)isnumeric(x),varargsout);
    if ~all(numericOutput)
        error('Non-numeric output arguments by the input function are not accepted.');
    end
    
    % Loop over the different outputs of the function
    for iOut = 1:nout
        out = varargsout{iOut};
        
        if isempty(sizeOut{iOut})
            sizeOut{iOut} = size(out);
        elseif any(sizeOut{iOut}~=size(out))
            error(['Inconsistent output variable size. ',...
                'One of the outputs of the analyzed function is changing its size in between runs. ',...
                'To solve this, fix the axis of the output and interpolate the result. \n%s'],...
                '  Ex/ outFix = interp1(varAxis,out,fixAxis,''pchip'')')
        end
        % Convert vectors to rows
        if iscolumn(out)
            out = out.';
        end
        nParam(iOut) = numel(out);
        % Store outputs in an N-dimensional array
        evals{iOut} = cat(1,evals{iOut},shiftdim(out,-1));
    end
    
    if Verbose
        fprintf(repmat('\b',1,numel(S)));
    end
end
warning('on','all')

% Compile statistics for all parameters from bootstrap samples
%-------------------------------------------------------------------------------
for iOut = 1:numel(evals)
    
    boots = evals{iOut};
    means = mean(boots);
    medians = median(boots);
    stds = std(boots);
    
    for i = 1:nParam(iOut)
        booti = boots(:,i);
        
        %Statistical metrics
        stats{iOut}(i).mean = means(i);
        stats{iOut}(i).median = medians(i);
        stats{iOut}(i).std = stds(i);
        stats{iOut}(i).p2 = percentile(booti,2,1);
        stats{iOut}(i).p25 = percentile(booti,25,1);
        stats{iOut}(i).p75 = percentile(booti,75,1);
        stats{iOut}(i).p98 = percentile(booti,98,1);
        if nParam(iOut)<20
            %Kernel-density estimations
            xmin = 0.9*stats{iOut}(i).p2;
            xmax = 1.1*stats{iOut}(i).p98;
            if all(diff(booti)==0)
                y = 1;
                x = 0;
            else
                [~,y,x] = kde(booti,100,xmin,xmax);
            end
            stats{iOut}(i).bootdist.x = x;
            stats{iOut}(i).bootdist.y = y;
            edges = linspace(xmin,xmax,50);
            
            %Histograms
            %Determine optimal bins via Freedman-Diaconis rule
            nbins = round(range(booti)/(2*iqr(booti)/(numel(booti)).^(1/3)),0);
            if isinf(nbins) || isnan(nbins)
                nbins = 10;
            end
            
            [bins,edges] = histcounts(booti,nbins,'Normalization','pdf');
            stats{iOut}(i).boothist.bins = bins;
            stats{iOut}(i).boothist.edges = edges;
        end
    end
    
end

end
