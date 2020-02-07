%
% PREPVALIDATION Computes parameter combinations for sensitivity analysis
% 
%   varparam = PREPVALIDATION(param)
%   Returns all the possible permutations of the parameters in the input
%   structure array (param). All possible combinations are randomly permuted and 
%   returned as a cell array. The (param) structure must have fields
%   corresponding to the names of the parameters. The fields contain all
%   the values/strings/logicals which the parameter can adopt.
% 
%   varparam = PREPVALIDATION(param,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
%       'randperm' - Specifies whether to randomly permute the combination
%                    vectors of validation parameter values (default = true)
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.


function varparam = prepvalidation(Parameters,varargin)

if ~isa(Parameters,'struct')
    error('First input must be a structure.');
end

% Parse and validate optional input
randomizeOrdering = parseoptional({'randperm'},varargin);

if isempty(randomizeOrdering)
    randomizeOrdering = true;
else
    validateattributes(randomizeOrdering,{'logical'},{'nonempty'},mfilename,'randperm');
end

% Get parameter names, values, etc
ParNames = fieldnames(Parameters);
ParValues = struct2cell(Parameters);
nParams = numel(ParNames);
nValues = cellfun(@numel,ParValues);
nCombinations = prod(nValues);

% Compile list of parameter combinations
[idx{1:nParams}] = ind2sub(nValues,1:nCombinations);
idx = cell2mat(idx.');

% Generate the parameter combinations
varparam = cell(nCombinations,nParams);
for c = 1:nCombinations
    for p = 1:nParams
        FieldValues = ParValues{p};
        if iscell(FieldValues)
            varparam{c,p} = FieldValues{idx(p,c)};
        else
            varparam{c,p} = FieldValues(idx(p,c));
        end
    end
end

% Randomize the order of the list
if randomizeOrdering
    varparam = varparam(randperm(nCombinations),:);
end
