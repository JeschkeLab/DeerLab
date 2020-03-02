function [pass,maxerr] = test(opt)

% Check that metric distances are computed allright

rng(1);
trial = rand(20,1);
trial = trial/sum(trial);
truth = trial;
metricNames = {'overlap','determination','chebyshev','cosine',...
    'correlation','chi','bregman','mad','msd','rmsd', 'nrmsd',...
    'hellinger','euclidean', 'braycurtis','bhattacharyya','tv'};
for i = 1:length(metricNames)
    values(i) = metrics(trial,truth,metricNames{i});
end
values_all = metrics(trial,truth);

% Pass 1: all metrics identify that values are equal
pass(1)  = all(values < 1e-7);
% Pass 1: all metrics are real values
pass(2)  = all(isreal(values));
% Pass 1: full metric is returned as structure
pass(3) = isstruct(values_all);

pass = all(pass);
maxerr = NaN;
 

end

