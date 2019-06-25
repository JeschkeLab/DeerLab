function [S,scale] = deer_sim(r,distr,t,exci)
% Simulation of DEER signal sim(t) from distance distribution distr(r)
% and scaling factor sc of distribution so that it fits original data set dipevo(t)
%
% Input:
%    r      distance vector, in nm
%    distr  distance distribution
%    t      time vector, in us
%    exci   excitation bandwidth, in MHz
% Output:
%    S      time-domain DEER signal
%    scale  scaling factor

bwlimit = exist('exci','var');

ny0 = 52.04; % dipolar frequency at 1 nm (MHz)
threshold = 1e-3*max(distr); % minimimum level for simulation of a point

S = zeros(size(t));
mdepth = 0;
fdepth = 0;
neglect = distr<=threshold;
if ~bwlimit, weight = 1; end
nTheta = 1000;
costheta = (0:nTheta-1)/nTheta; % cos(theta)
for k = 1:length(r)
    % simulate only for distances with significant contribution
    if neglect(k), continue; end
    
    nydd = ny0/r(k)^3; % dipolar frequency at current distance
    nyac = nydd*(3*costheta.^2-1); % dipolar frequency for all orientations
    for m = 1:nTheta % theta loop for orientation dependence
        if bwlimit
            weight = exp(-(nyac(m)/exci)^2);
            mdepth = mdepth + weight*distr(k);
        end
        fdepth = fdepth + distr(k);
        S = S + distr(k)*cos(2*pi*nyac(m)*t)*weight; % add contribution to deer signal
    end
end
if bwlimit
    scale = mdepth/fdepth;
    S = 0.01*S/mdepth;
else
    scale = 1;
    S = 0.01*S/fdepth;
end
S = S + 0.99*ones(size(S));
