% Runs DEER processing using the specified network file. Syntax:
%
%       answers=process_using(deer_traces,net_file_name)
%
% Parameters:
%
%   deer_traces  -  DEER trace(s) as a column vector or a 
%                   matrix with multiple columns. The num-
%                   ber rows must match the input size of
%                   the neural network.
%
%   net_file_name - the name of the ceckpoint file produ-
%                   ced by the train() function of Matlab
%                   Neural Network Toolbox, including the
%                   full path and the .mat extension.
%
% Outputs:
%
%   answers       - distance distribution(s) as a column
%                   vector or a matrix with multiple co-
%                   lumns. The number of rows matches the
%                   output size of the neural network.
%
% Note: the function assumes a simple feed-forward topology and
%       will fail if that is not the case.
%
% i.kuprov@soton.ac.uk
% s.g.worswick@soton.ac.uk

function answers=process_using(deer_traces,net_file_name)

% Check consistency
grumble(deer_traces,net_file_name);

% Try to load the network
load(net_file_name,'checkpoint');
if ~exist('checkpoint','var')
    error([net_file_name ' does not contain a neural network checkpoint.']);
end

% Check that input dimensions match
if size(checkpoint.net.IW{1},2)~=size(deer_traces,1)
    error([net_file_name 'input layer size does not match data size.']);
end

% % Define transfer functions if Deep Learning Toolbox is not installed
% if ~exist('tansig','file') % tansig is from Deep Learning Toolbox
%   tansig = @(n) 2./(1+exp(-2*n))-1;
% end
% if ~exist('logsig','file') % logsig is from Deep Learning Toolbox
%   logsig = @(n) 1./(1+exp(-n));
% end

% Apply the input layer
answers=tansig(checkpoint.net.IW{1}*deer_traces+checkpoint.net.b{1});
    
% Apply middle layers
for n=2:(checkpoint.net.numLayers-1)
    answers=tansig(checkpoint.net.LW{n,n-1}*answers+checkpoint.net.b{n}); 
end

% Apply the last layer
if strcmp(checkpoint.net.layers{end}.transferFcn,'logsig')
    answers=logsig(checkpoint.net.LW{end,end-1}*answers+checkpoint.net.b{end}); 
elseif strcmp(checkpoint.net.layers{end}.transferFcn,'tansig')
    answers=tansig(checkpoint.net.LW{end,end-1}*answers+checkpoint.net.b{end});
else
    error('unexpected transfer function at the last layer.');
end
    
end

% Consistency enforcement
function grumble(deer_traces,net_file_name)
if (~isnumeric(deer_traces))||(~isreal(deer_traces))||...
   (~ismatrix(deer_traces))
    error('deer_traces must be a real column vector or matrix.');
end
if (~ischar(net_file_name))||(~strcmp(net_file_name((end-3):end),'.mat'))
    error('net_file_name must be a character string ending in .mat');
end
end

% Those who say that all cultures are equal never explain why 
% the results of those cultures are so grossly unequal.
%
% Thomas Sowell

