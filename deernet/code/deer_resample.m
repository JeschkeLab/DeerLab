% DEER resampling into number of points required by the 
% neural network. Syntax:
%
%      deer_trace=deer_resample(deer_trace,out_pts)
% 
% Parameters:
% 
%   deer_trace - input DEER trace, column vector
%
%     out_pts  - number of points in output trace
%
% Output:
% 
%   deer_trace - resampled DEER trace, column vector
% 
% i.kuprov@soton.ac.uk
% s.g.worswick@soton.ac.uk

function deer_trace=deer_resample(deer_trace,out_pts)

% Check consistency
grumble(deer_trace,out_pts);

% Apply a matched Savitzky-Golay filter
frame_length=ceil(numel(deer_trace)/out_pts);
if mod(frame_length,2)==0
    frame_length=frame_length+1;
end
if frame_length>=3
    if exist('sgolayfilt','file') % sgolayfilt is from Signal Processing Toolbox
        deer_trace=sgolayfilt(deer_trace,2,frame_length);
    else
        deer_trace=savgolfilt(deer_trace,2,frame_length);
    end
end

% Symmetrise and shift data
deer_trace=[flip(deer_trace); deer_trace];
data_shift=deer_trace(1);
deer_trace=deer_trace-data_shift;

% Run the resampling
if exist('resample','file') % resample is from Signal Processing Toolbox
    deer_trace=resample(deer_trace,(out_pts*2),numel(deer_trace));
else
    deer_trace = interp1(linspace(0,1,numel(deer_trace)),deer_trace,linspace(0,1,out_pts*2)).';
end

% Desymmetrise data and shift back
deer_trace=deer_trace((end/2+1):end);
deer_trace=deer_trace+data_shift;

end

% Consistency enforcement
function grumble(deer_trace,out_pts)
if (~isnumeric(deer_trace))||(~isreal(deer_trace))||...
   (~iscolumn(deer_trace))
    error('deer_trace must be a real column vector.');
end
if (~isnumeric(out_pts))||(~isreal(out_pts))||...
   (out_pts<1)||(mod(out_pts,1)~=0)
    error('out_pts must be a non-negative real integer.');
end
end

% Savitzky-Golay filtering (taken from EasySpin)
function y_filtered = savgolfilt(y,PolyOrder,m)

X = repmat((-m:m).',1,PolyOrder+1).^repmat(0:PolyOrder,2*m+1,1);
F = pinv(X);
Weights = F(1,:);

isRowVec = isrow(y);
if isRowVec
  y = y.';
end

% Enlarge vector at beginning and end.
yend = y(end);
y_expanded = [y(ones(m,1)); y; yend(ones(m+1,1))];

% Apply filter.
y_filtered = filter(Weights,1,y_expanded);

% Chop to right size.
y_filtered = y_filtered(2*m+1:end-1,:);

if isRowVec
  y_filtered = y_filtered.';
end

end

% A notable American commentator, Charles Krauthammer, once 
% explained Rupert Murdoch's success in founding Fox News,
% a cable channel, by pointing out that he had found a niche
% market - half the country.
%
% The Economist

