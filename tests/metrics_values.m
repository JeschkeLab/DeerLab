function [err,data] = test(opt,olddata)

%======================================================
% Make sure metric distances are computed allright
%======================================================

metricNames = {'overlap','determination','chebyshev','cosine',...
    'correlation','chi','bregman','mad','msd','rmsd', 'nrmsd',...
    'hellinger','euclidean', 'bhattacharyya','tv'};

rng(2)
trial = rand(20,1);
trial  = trial/sum(trial);
truth = trial;

for i=1:length(metricNames)
values(i) = metrics(trial,truth,metricNames{i});
end

err(1)  = any(values>1e-7);
err(2)  = any(~isreal(values));
err = any(err);

data = [];


end

