%
% PREPVALIDATION Computes the validation parameter combinations
% 
%   varparam = PREPVALIDATION(param)
%   Returns all the possible permutations of the parameters in the input
%   structure array (param). All possible combinations are randomly permuted and 
%   returned as a cell array. The (param) structure must have to following
%   fields for all of the N parameters:
%
%           param(N).name - Name of the parameter variable
%           param(N).values - Values or strings to be evaluated (array or
%                             cell array)
% 
%   varparam = PREPVALIDATION(param,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
%       'randperm' - Specifies whether to randomly permute the combination
%                    vectors of validation parameter values (default = true)
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.

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
