function [err,data] = test(opt,olddata)

%======================================================
% Make sure metric distances are computed allright
%======================================================

metricNames = {'overlap','determination','chebyshev','cosine',...
    'correlation','chi','bregman','mad','msd','rmsd', 'nrmsd',...
    'hellinger','euclidean', 'bhattacharyya','tv'};

trial = rand(20,1);
trial  = trial/sum(trial);
truth = trial;

for i=1:length(metricNames)
values(i) = metrics(trial,truth,metricNames{i});
end

err  = any(values>1e-10);
err  = any(~isreal(values));

data = [];


end

