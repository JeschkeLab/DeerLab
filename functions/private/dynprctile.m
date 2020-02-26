%
% DYNPRCTILE Dynamic percentile estimation
%
%    [perc,markers] = DYNPRCTILE(p,sample,markers)
%    Estimates the pth-percentile of a dynamic set of samples. For the
%    first call the markers structure must be passed empty, and it will be
%    returned filled. For each additional sample added, the filled markers
%    structure must be passed back.
%
%    This function is a multi-dimensional generalization of the original
%    P2-algorithm by Jain and Chlamtac (see function below). 

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function [perc,markers] = dynprctile(p,sample,markers)

%Store the dimensionality of the input sample
origSize = size(sample);

%Fold into single column vector
sample = sample(:);

%Get number of points in sample
N = numel(sample);

%Check if the markers structure has already been cosntructed
if isempty(markers)
    n = cell(1,N);
    np = cell(1,N);
    n(:) = {1:5};
    np(:) = {[1, 1 + 2*p, 1+4*p, 3 + 2*p, 5]};
    q = cell(1,N);
else
    n = markers.n;
    np = markers.np;
    q = markers.q;
end

%For the first 5 samples, just fill the markers
if numel(q{1})<5
    for j=1:numel(sample)
        q{j}(end+1) = sample(j);
    end
    perc = sample;
else
    %For the rest of samples, run the P2-algorithm to update the percentile
    idx = 1:N;
    [perc,n,np,q] = arrayfun(@(idx) p2algorithm(p,sample(idx),n{idx},np{idx},q{idx}),idx,'UniformOutput',false);
    %Convert back cell array returned by arrayfun
    perc = cell2mat(perc);
end
%Store status of markers for next function call
markers.n = n;
markers.np = np;
markers.q = q;

%Reshape the percentile back to the original sample dimensionality
perc = reshape(perc,origSize);

%Ensure columns for the vectors
if isvector(perc) 
    perc = perc(:);
end

end

function [quant,n,np,q] = p2algorithm(p,sample,n,np,q)
%
% P2--Algorithm for Dynamic Calculation of Quantiles
%
% See: 
%
% Raj Jain and Imrich Chlamtac. 1985. 
% The P2 algorithm for dynamic calculation of quantiles and histograms 
% without storing observations. 
% Commun. ACM 28, 10 (October 1985), 1076–1085. 
% DOI: 
%     https://doi.org/10.1145/4372.4378
%
% Implemented by Luis Fabregas, 2020
% NOTES: There are some typos in the original paper
%
%  - In the Parabolic Prediction formula, the last line (".n1 <- ni + d")
%    there is no multiplication, just a new line. The actual formula ends
%    on the second line.
%   
%  - In Box 1, the section B.2. the increment should say: "ni <- ni+1 i=k+1,...,5"
%
%  - In Box 1, the section B.3. the in the IF statement it should read
%    "IF qi-1 < qiprime < qi+1"
%    "THEN qi <- qiprime"
%

dnp = [0, p/2, p, (1+p)/2, 1];

%Find cell for current sample, and adjust extreme values
x = sample;
if x < q(1)
    k = 1;
    q(1) = x;
elseif q(1)<= x && x < q(2)
    k = 1;
elseif q(2)<= x && x < q(3)
    k = 2;
elseif q(3)<= x && x < q(4)
    k = 3;
elseif q(4)<= x && x <= q(5)
    k = 4;
elseif q(5)< x
    k = 4;
    q(5) = x;
end

%Increment marker positions
n(k+1:5) = n(k+1:5) + 1;
np = np + dnp;

%Adjust heights of middle markers
for i = 2:4
    d = np(i) - n(i);
    if (d >= 1 && (n(i+1) - n(i)) > 1) || (d <= -1 && (n(i-1) - n(i)) < -1)
        d = sign(d);
        %Try adjusting height using parabolic prediction (P2) formula
        qp = q(i) + d/(n(i+1) - n(i-1))*( (n(i) - n(i-1) + d)*(q(i+1) - q(i))/(n(i+1) - n(i)) + (n(i+1) - n(i) - d)*(q(i) - q(i-1))/(n(i) - n(i-1)) );
        if q(i-1)<qp && qp<q(i+1)
            q(i) = qp;
        else
            q(i) = q(i) + d*(q(i+d) - q(i))/(n(i+d) - n(i));
        end
        n(i) = n(i) + d;
    end
end

%Return the requested quantile
quant = q(3);

end