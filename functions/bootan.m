function stats = bootan(fcn,V,Vfit,varargin)


% Parse the optional parameters in varargin
optionalProperties = {'Verbose','Samples'};
[Verbose,Nsamples] = parseoptional(optionalProperties,varargin);


if isempty(Nsamples)
    Nsamples = 1000;
else
    validateattributes(Nsamples,{'numeric'},{'nonnegative','scalar','nonempty'},mfilename,'Nsamples')
end

if isempty(Verbose)
    Verbose = false;
else
    validateattributes(Verbose,{'logical'},{'nonempty'},mfilename,'Verbose')
end


% Determine the number of outputs returned by function handle and allocate arrays
nout = getmaxnargout(fcn,V);
evals = cell(1,nout);
sizeOut = cell(1,nout);

%Get residual
res = Vfit - V;
N = numel(res);

warning('off','all')
for i=1:Nsamples
    
    %Inform of progress if requested
    if Verbose
        S = sprintf('Bootstrapping: %i/%i samples finished',i,Nsamples);
        fprintf(S);
    end
    
    if i>1
        % Get a new bootstrap sample
        sampling = randperm(N);
        Vboot = Vfit + res(sampling);
    else
        Vboot = V;
    end
    
    % Run the user function with current factor set
    varargsout = cell(1,nout);
    [varargsout{:}] = fcn(Vboot);
    
    % Assert that all outputs are strictly numerical
    numericOutput = cellfun(@(x)isnumeric(x),varargsout);
    if ~all(numericOutput)
        error('Non-numeric output arguments by the input function are not accepted.');
    end
    
    %Loop over the different outputs of the function
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
        nParam(iOut) = numel(out);
        % Store outputs in an N-dimensional array
        evals{iOut} = cat(1,evals{iOut},shiftdim(out,-1));
    end
    
    if Verbose && exist('S','var')
        fprintf(repmat('\b',1,numel(S)));
    end
end
warning('on','all')

% Loop over all fitted parameters

for iOut = 1:numel(evals)
    boots = evals{iOut};
for i=1:nParam(iOut)
    
    stats{iOut}(i).mean = mean(boots(:,i));
    stats{iOut}(i).median = median(boots(:,i));
    stats{iOut}(i).std = std(boots(:,i));
    [~,y,x] = kde(boots(:,i),400); 
    stats{iOut}(i).bootdist.x = x; 
    stats{iOut}(i).bootdist.y = y;
    stats{iOut}(i).p2 = percentile(boots(:,i),2,1);
    stats{iOut}(i).p25 = percentile(boots(:,i),25,1);
    stats{iOut}(i).p75 = percentile(boots(:,i),75,1);
    stats{iOut}(i).p98 = percentile(boots(:,i),98,1);
end

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