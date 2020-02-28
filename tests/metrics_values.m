function [pass,maxerr] = test(opt)

%======================================================
% Make sure metric distances are computed allright
%======================================================

metricNames = {'overlap','determination','chebyshev','cosine',...
    'correlation','chi','bregman','mad','msd','rmsd', 'nrmsd',...
    'hellinger','euclidean', 'braycurtis','bhattacharyya','tv'};

rng(2);

trial = rand(20,1);
trial = trial/sum(trial);
truth = trial;

for i = 1:length(metricNames)
    values(i) = metrics(trial,truth,metricNames{i});
end
values_all = metrics(trial,truth);

err(1)  = any(values>1e-7);
err(2)  = any(~isreal(values));
err(3) = ~isstruct(values_all);

pass = all(err);
maxerr = NaN;
 

end

