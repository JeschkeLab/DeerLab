function metric = metrics(pointA,pointB,metricName) 
if nargin<3
   metricName = 'all'; 
end

  switch metricName
      case 'overlap'
          metric = 1 - sum(min(pointA,pointB));
      case 'determination'
          metric = sum((pointA - pointB).^2)/sum((pointA - mean(pointA)).^2);
      case 'chebyshev'
          metric = max(abs(pointA - pointB));
      case 'braycurtis'
          metric = 1 - dot(pointA - mean(pointA),pointB - mean(pointB))/(norm(pointB - mean(pointB))*norm(pointA - mean(pointA)));
      case 'cosine'
          metric = norm(pointA - pointB,1)/norm(pointA + pointB,1);
      case 'correlation'
          metric = 1 - dot(pointA,pointB)/(norm(pointB)*norm(pointA));
      case 'chi'
          metric = sum((pointA - pointB).^2./(pointA + pointB));
      case 'bregman'
          metric = 1/2*norm(pointA - pointB) + dot(pointA, pointA - pointB);
      case 'mad'
          metric = mean(abs(pointA - pointB));
      case 'msd'
          metric = mean(abs(pointA - pointB).^2);
      case 'rmsd'
          metric = sqrt(mean(abs(pointA - pointB)));
      case 'nrmsd'
          metric = sqrt(mean(abs(pointA - pointB)))/(max(pointA) - min(pointA));
      case 'hellinger'
          metric = sqrt(1 - sum(sqrt(pointA.*pointB) / sqrt(sum(pointA) * sum(pointB))));
      case 'euclidean'
          metric  = norm(pointA - pointB);
      case 'bhattacharyya'
          metric  = -log(sum(sqrt(pointA.*pointB))/sqrt(sum(pointA))/sqrt(sum(pointB)));
      case 'tv'
          metric =  1/2*norm(pointA - pointB,1);
      otherwise
          metric.Overlap = sum(min(pointA/sum(pointA),pointB/sum(pointB)));
          metric.Determination = sum((pointA - pointB).^2)/sum((pointA - mean(pointA)).^2);
          metric.Chebyshev = max(abs(pointA - pointB));
          metric.BrayCurtis = norm(pointA - pointB,1)/norm(pointA + pointB,1);
          metric.Cosine = 1 - dot(pointA,pointB)/(norm(pointB)*norm(pointA));
          metric.Correlation = 1 - dot(pointA - mean(pointA),pointB - mean(pointB))/(norm(pointB - mean(pointB))*norm(pointA - mean(pointA)));
          metric.Chi = sum((pointA - pointB).^2./(pointA + pointB));
          metric.Bregman = 1/2*norm(pointA - pointB) + dot(pointA, pointA - pointB);
          metric.MeanAbsDev = mean(abs(pointA - pointB));
          metric.MeanSquaredDev = mean(abs(pointA - pointB).^2);
          metric.RMSDdeviation = sqrt(mean(abs(pointA - pointB)));
          metric.normalizedRMSDdeviation = sqrt(mean(abs(pointA - pointB)))/(max(pointA) - min(pointA));
          metric.Euclidean  = norm(pointA - pointB);
          metric.Hellinger = sqrt(1 - sum(sqrt(pointA.*pointB) / sqrt(sum(pointA) * sum(pointB))));
          metric.TotalVariation =  1/2*norm(pointA - pointB,1);
          metric.Bhattacharyya  = -log(sum(sqrt(pointA.*pointB))/sqrt(sum(pointA))/sqrt(sum(pointB)));         
  end
end