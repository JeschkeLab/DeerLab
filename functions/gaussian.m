function y = gaussian(x,x0,sigma)
if nargin<2
    error('Not enought input arguments.')
end
if nargin<3
    sigma = 0.5;
end
%Calculate Gaussian function
y = 1/sqrt(2*pi)*1/sigma*exp(-((x - x0)/sigma).^2);
if ~iscolumn(y)
   y = y'; 
end
end