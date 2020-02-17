function metric = metrics(A,B,metricName) 

if nargin<3
   metricName = 'all'; 
end

switch metricName
    case 'overlap'
        metric = 1 - sum(min(A,B));
    case 'determination'
        metric = sum((A - B).^2)/sum((A - mean(A)).^2);
    case 'chebyshev'
        metric = max(abs(A - B));
    case 'braycurtis'
        metric = 1 - dot(A - mean(A),B - mean(B))/(norm(B - mean(B))*norm(A - mean(A)));
    case 'cosine'
        metric = norm(A - B,1)/norm(A + B,1);
    case 'correlation'
        metric = 1 - dot(A,B)/(norm(B)*norm(A));
    case 'chi'
        metric = sum((A - B).^2./(A + B));
    case 'bregman'
        metric = 1/2*norm(A - B) + dot(A, A - B);
    case 'mad'
        metric = mean(abs(A - B));
    case 'msd'
        metric = mean(abs(A - B).^2);
    case 'rmsd'
        metric = sqrt(1/numel(A)*norm(A - B)^2);
    case 'nrmsd'
        metric = sqrt(mean(abs(A - B)))/(max(A) - min(A));
    case 'hellinger'
        metric = sqrt(1 - sum(sqrt(A.*B) / sqrt(sum(A) * sum(B))));
    case 'euclidean'
        metric  = norm(A - B);
    case 'bhattacharyya'
        metric  = -log(sum(sqrt(A.*B))/sqrt(sum(A))/sqrt(sum(B)));
    case 'tv'
        metric =  1/2*norm(A - B,1);
    otherwise
        metric.Overlap = metrics(A,B,'overlap'); 
        metric.Determination = metrics(A,B,'determination');
        metric.Chebyshev = metrics(A,B,'chebyshev'); 
        metric.BrayCurtis = metrics(A,B,'braycurtis');
        metric.Cosine = metrics(A,B,'cosine');
        metric.Correlation = metrics(A,B,'correlation');
        metric.Chi = metrics(A,B,'chi');
        metric.Bregman = metrics(A,B,'bregman');
        metric.MeanAbsDev = metrics(A,B,'mad');
        metric.MeanSquaredDev = metrics(A,B,'msd');
        metric.RMSDdeviation = metrics(A,B,'rmsd');
        metric.normalizedRMSDdeviation = metrics(A,B,'nrmsd');
        metric.Euclidean = metrics(A,B,'euclidean');
        metric.Hellinger = metrics(A,B,'hellinger');
        metric.Bhattacharyya  = metrics(A,B,'bhattacharyya');
        metric.TotalVariation =  metrics(A,B,'tv');
end

end
