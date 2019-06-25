% The six examples used in the paper by Worswick, Spencer, 
% Jeschke, and Kuprov (in press).
%
% s.g.worswick@soton.ac.uk
% gunnar.jeschke@phys.chem.ethz.ch
% i.kuprov@soton.ac.uk

function example_set_paper()

% Sample I, networks optimized for sharp peaks
netset_location='net_set_sharp_peaks';
expt_data=load('sample_set_paper/sample_I_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),netset_location); drawnow();

% Sample II, networks optimized for broad peaks
netset_location='net_set_broad_peaks';
expt_data=load('sample_set_paper/sample_II_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),netset_location); drawnow();

% Sample III, networks optimized for sharp peaks
netset_location='net_set_sharp_peaks';
expt_data=load('sample_set_paper/sample_III_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),netset_location); drawnow();

% Sample IV, networks optimized for broad peaks
netset_location='net_set_broad_peaks';
expt_data=load('sample_set_paper/sample_IV_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),netset_location); drawnow();

% Sample V, networks trained for all peak types
netset_location='net_set_any_peaks';
expt_data=load('sample_set_paper/sample_V_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),netset_location); drawnow();

% Sample VI, networks optimized for broad peaks
netset_location='net_set_broad_peaks';
expt_data=load('sample_set_paper/sample_VI_DEERNet_input.dat','-ASCII');
deernet(expt_data(:,2),1e-6*expt_data(:,1),netset_location); drawnow();

end

