function varparam = prepvalidation(Parameters,varargin)

%Parse & validate required input
if ~isfield(Parameters,{'name','values'})
   error('The input structure must contain the ''name'' and ''values'' fields.') 
end

%Parse & validate optional input
[randpermflag] = parseoptional({'randperm'},varargin);

if isempty(randpermflag)
    randpermflag = true;
else
    validateattributes(randpermflag,{'logical'},{'nonempty'},mfilename,'randperm')
end


%Get number of variables to validate
nParam = numel(Parameters);
%Get number of trials for each validation parameter
varTrials = zeros(nParam,1);
for i=1:nParam
varTrials(i) = length(Parameters(i).values);
end
totalTrials = prod(varTrials);
%Preallocate validation parameters
varparam = cell(totalTrials,nParam);
%Loop over all validation trials
for Pindex = 1:totalTrials
    idx=cell(1,numel(varTrials));
    [idx{:}] = ind2sub(varTrials,Pindex);
    idx = cell2mat(idx);
    %Generate the validation parameter values
    for varIdx = 1:nParam
        %If user supplies the vector directly then take it from there
        if iscell(Parameters(varIdx).values)
            value = Parameters(varIdx).values{idx(varIdx)};
        else
            value = Parameters(varIdx).values(idx(varIdx));
        end
        %Save current sampled value
        varparam{Pindex,varIdx} = value;
    end
end

%Shuffle the order of parameter combination sets
if randpermflag
    varparam = varparam(randperm(size(varparam,1)),:);
end
