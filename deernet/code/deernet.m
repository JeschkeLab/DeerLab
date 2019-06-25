% Use pre-trained neural networks to extract distance distributions from
% primary DEER data for regular cases. Syntax:
% 
%           [dist_axis,outcomes]=deernet(deer_trace,time_axis,net_dir,silent)
% 
% Parameters:
% 
%    deer_trace - experimental DEER trace, a column vector
%
%    time_axis  - experimental time axis, a column vector in 
%                 microseconds, must start at zero and have
%                 a uniform time step
%
%    net_dir    - directory, where the networks are stored    
%
%    silent     - supresses output to command window if present and true
%                 Boolean flag
%
%    net_dir    - directory containing the network ensemble
%                 as multiple *.m files produced by genFunction
%                 routine of Matlab's Neural Network toolbox.
%                 Network file names must start with "n_".
%
%    parameters - training parameters used to generate the
%                 network ensemble
% 
% Output:
% 
%    dist_axis  - distance axis, a column vector (Angstrom)
%
%    outcomes   - distance distributions, a horizontal array
%                 of column vectors. The number of column
%                 vectors is equal to the number of networks
%                 in the ensemble.
% 
% s.g.worswick@soton.ac.uk
% gunnar.jeschke@phys.chem.ethz.ch
% i.kuprov@soton.ac.uk

function [dist_axis,outcomes]=deernet(deer_trace,time_axis,net_dir,silent)

% Check consistency
grumble(deer_trace,time_axis,net_dir);

% Get the parameter set
if ~exist([net_dir filesep 'netset_params.m'],'file')
    error('netset_params.m file is missing from the directory specified.');
else
    run([net_dir filesep 'netset_params.m']);
end

if ~exist('silent','var') || isempty(silent)
    silent = false;
end

% Scale and shift the input data
deer_trace=deer_trace-deer_trace(end);
deer_trace=deer_trace/deer_trace(1);

% Resample the input data
deer_trace=deer_resample(deer_trace,parameters.npoints);
time_axis=linspace(time_axis(1),time_axis(end),parameters.npoints);

% Read the ensemble
if ~exist([net_dir filesep 'good_nets.mat'],'file')
    error('good_nets.mat file is missing from the directory specified.');
else
    load([net_dir filesep 'good_nets.mat'],'good_nets');
    good_nets=good_nets; %#ok<ASGSL,NODEF>
end

% Run the network ensemble
outcomes=zeros(size(deer_trace,1),numel(good_nets));
parfor n=1:numel(good_nets)
    
    % Update the user
    disp(['Network ' num2str(n) ' out of ' num2str(numel(good_nets)) '...']);
    
    % Run each network
    outcomes(:,n)=process_using(deer_trace,[net_dir filesep good_nets{n}]);
    
end

% Rescale the output distribution for experimental tmax
scale=(time_axis(end)/parameters.max_time)^(1/3);
dist_axis=linspace(parameters.min_dist*scale,...
                   parameters.max_dist*scale,parameters.npoints)';
               
% If no outputs requested, do the plotting
if nargout==0
    
    % Plot the experimental data
    figure(); subplot(1,3,1); plot(1e6*time_axis,deer_trace);
    axis tight; title('experimental data'); grid on;
    xlabel('time, $\mu$s','interpreter','LaTex');
    ylabel('amplitude, a.u.','interpreter','LaTex');
    
    % Plot the distance distribution ensemble
    subplot(1,3,2); plot(dist_axis,outcomes);
    title('network ensemble'); axis tight; grid on;
    xlabel('distance, $\rm{\AA}$','interpreter','LaTex');
    
    % Plot the average distance distribution
    mean_outcome=mean(outcomes,2);
    subplot(1,3,3); hold on;
    plot(dist_axis,mean_outcome);
    
    % Run the statistics
    lower_bound=mean_outcome-2*std(outcomes,[],2);
    lower_bound(lower_bound<0)=0;
    upper_bound=mean_outcome+2*std(outcomes,[],2);
    
    % Plot the statistics
    area(dist_axis,upper_bound,'FaceAlpha',0.25,'EdgeAlpha',0);
    area(dist_axis,lower_bound,'FaceColor',[1 1 1],'EdgeAlpha',0);
    axis tight; grid on; title('ensemble statistics');
    set(gca,'Layer','top'); legend({'mean','95%'});
    xlabel('distance, $\rm{\AA}$','interpreter','LaTex');
    pos=get(gcf,'pos'); set(gcf,'pos',[pos(1) pos(2) 900 300]);
    
end

end

% Consistency enforcement
function grumble(deer_trace,time_axis,net_dir)
if (~isnumeric(deer_trace))||(~isreal(deer_trace))||...
   (~iscolumn(deer_trace))
    error('deer_trace must be a real column vector.');
end
if (~isnumeric(time_axis))||(~isreal(time_axis))||...
   (~iscolumn(time_axis))
    error('time_axis must be a real column vector.');
end
if numel(deer_trace)~=numel(time_axis)
    error('deer_trace and time_axis must have the same length.');
end
if ~ischar(net_dir)
    error('net_dir must be a character string.');
end
end

% When will colleges start to take intellectual and political
% diversity as seriously as they take the more superficial 
% forms of diversity? 
%
% Milo Yiannopoulos

