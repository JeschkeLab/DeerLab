% Processing of DEER traces in the presence of significant exchange
% interaction. Data from Figure 10 in the Supplementary Information
% of http://dx.doi.org/10.1038/ncomms14842, very kindly provided by
% Sabine Richert.
%
% s.g.worswick@soton.ac.uk
% gunnar.jeschke@phys.chem.ethz.ch
% richert.sabine@gmail.com
% i.kuprov@soton.ac.uk

function example_set_exchange()

% Exchange-resistant networks
netset_location='net_set_with_exch';

% Exchanging system - doubly linked
expt_data=load('sample_set_exchange/fig10_SI_P2P2_time_phasecorrDEERxy.txt','-ASCII');
positive_times=(expt_data(:,1)>=0); expt_data=expt_data(positive_times,:);
deernet(expt_data(1:200,2),1e-6*expt_data(1:200,1),netset_location); drawnow();

% Exchanging system - singly linked
expt_data=load('sample_set_exchange/fig10_SI_P2X_time_phasecorrDEERxy.txt','-ASCII');
positive_times=(expt_data(:,1)>=0); expt_data=expt_data(positive_times,:);
deernet(expt_data(1:200,2),1e-6*expt_data(1:200,1),netset_location); drawnow();

% Non-exchanging system - unlinked
expt_data=load('sample_set_exchange/fig10_SI_XX_time_phasecorrDEERxy.txt','-ASCII');
positive_times=(expt_data(:,1)>=0); expt_data=expt_data(positive_times,:);
deernet(expt_data(:,2),1e-6*expt_data(:,1),netset_location); drawnow();

% Standard networks
netset_location='net_set_any_peaks';

% Non-exchanging system - unlinked
expt_data=load('sample_set_exchange/fig10_SI_XX_time_phasecorrDEERxy.txt','-ASCII');
positive_times=(expt_data(:,1)>=0); expt_data=expt_data(positive_times,:);
deernet(expt_data(:,2),1e-6*expt_data(:,1),netset_location); drawnow();

end

