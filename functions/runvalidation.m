function runvalidation(ExpSignal,Parameters,varargin)

[Callback] = parseoptional({'Callback'},varargin);

if ~isempty(Callback)
    if ~isa(Callback,'function_handle')
        error('''Callback'' option must be a valid function handle.')
    end
end

varparam = prepareValidation(Parameters,'randperm',true);
varnames = {Parameters.name};
options = daoptions();

nParam = size(varparam,1);

Distributions = zeros(nParam,length(ExpSignal));

for i=1:nParam
    
    for j = 1:length(varnames) 
    options = setoptions(options,varnames{j},varparam{i,j});
    end
    %Run signal preparation
    ProcessedSignal = prepareSignal(ExpSignal,options);
    %Get the distribution
    Distributions(i,:) = getDistribution(ProcessedSignal,options);
    
    %Hook for GUI callbacks
    if ~isempty(Callback)
        Callback(i,Distributions);
    end
end

end
