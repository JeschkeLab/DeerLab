%
%   LEVENDIST Levenshtein distance between two strings
%
%   dist = LEVENDIST(str1,str2)
%   Computes the Levenshtein distance between strings (str1) and (str2)
%   using the Wagner–Fischer algorithm.
%
%   See: Wagner–Fischer algorithm
%        https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function dist = levendist(str1,str2)

%Get string length
n1 = 1 + numel(str1);
n2 = 1 + numel(str2);

%Prepare distance matrix
D = zeros(n1,n2);
D(:,1) = 0:n1-1;
D(1,:) = 0:n2-1;

%Flood filling of distance matrix
for i = 2:n1
  for j = 2:n2
      substcost = double(~(str1(i-1)==str2(j-1)));
      D(i,j) = min([D(i-1,j)+1, ... %deletion
                    D(i,j-1)+1, ... %insertion
                    D(i-1,j-1) + substcost]); %substitution
  end
end

%Extract Levenshtein distance
dist = D(end);

end