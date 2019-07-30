function varparam = prepareValidation(Parameters)

if ~isa(Parameters,'validationParameters')
    error('Input must be a valid validationParameters class object.')
end

%Get number of variables to validate
nParam = numel(Parameters);
%Get number of trials for each validation parameter
varTrials = [Parameters.trials];
totalTrials = prod([Parameters.trials]);
%Preallocate validation parameters
varparam = cell(totalTrials,nParam);
%Loop over all validation trials
for Pindex = 1:totalTrials
    idx=cell(1,numel(varTrials));
    [idx{:}] = ind2sub(varTrials,Pindex);
    idx = cell2mat(idx);
    %Generate the validation parameter values
    for varIdx = 1:nParam
        nPoints = Parameters(varIdx).trials;
        if isempty(Parameters(varIdx).vector)
            Range = Parameters(varIdx).range;
            %Sample the parameter values according to...
            if strcmp(Parameters(varIdx).sampling,'lin')
                %... linear sampling in the given range
                value =  Range(1) + idx(varIdx)*diff(Range)/nPoints;
            elseif strcmp(Parameters(varIdx).sampling,'rand')
                %... uniform random sampling in the given range
                randSamples = Range(1) +  diff(Range)*rand(1,nPoints);
                value = randSamples(idx(varIdx));
            end
        else
            %If user supplies the vector directly then take it from there
            if iscell(Parameters(varIdx).vector)
                value = Parameters(varIdx).vector{idx(varIdx)};
            else
                value = Parameters(varIdx).vector(idx(varIdx));
            end
        end
        %Save current sampled value
        varparam{Pindex,varIdx} = value;
    end
end