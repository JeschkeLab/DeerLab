%
%   ISNANOMETER
%
%   logical = ISNANOMETER(r)
%   If the distance axis is in nanometers, the function returns (true), if
%   it is in angstrom it returns (false).
%

function logical = isnanometer(r)

if iscolumn(r)
    r = r.';
end

if diff(minmax(r)) < 1 && length(r)>1
   warning('For the given distance axis (max(r) - min(r)) < 1. Assuming nanometer scale.') 
   logical = true; 
   return
end

if mean(r) < 2
   warning('For the given distance axis mean(r) < 2. Assuming nanometer scale.') 
   logical = true; 
   return
end

logical = log(mean(r)) < 2.99;


end