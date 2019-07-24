function varparam = prepareValidation(Parameters)

%Get number of variables to validate
nParam = numel(Parameters);
%Get number of trials for each validation parameter
varTrials = [Parameters.nPoints];
totalTrials = prod([Parameters.nPoints]);
%Preallocate validation parameters
varparam = zeros(totalTrials,nParam);
%Loop over all validation trials
for Pindex = 1:totalTrials
    idx=cell(1,numel(varTrials));
    [idx{:}] = ind2sub(varTrials,Pindex);
    idx = cell2mat(idx);
    %Generate the validation parameter values
    for varIdx = 1:nParam
        nPoints = Parameters(varIdx).nPoints;
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
        %Save current sampled value
        varparam(Pindex,varIdx) = value;
    end
end