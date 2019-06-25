% Processing of DEER traces in the presence of significant exchange
% interaction. Data from the example set of DeerAnalysis very kindly
% provided by Gunnar Jeschke.
%
% s.g.worswick@soton.ac.uk
% gunnar.jeschke@phys.chem.ethz.ch
% i.kuprov@soton.ac.uk

function example_set_deeran()

% Which net set directory to use
netset_location='net_set_any_peaks';

% Process the entire example set of DeerAnalysis
[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/CT_DEER_5nm/CT_DEER_5nm');
deernet(deer_trace,1e-6*time_axis,netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/CT_deer_broad/CT_deer_broad');
deernet(deer_trace,1e-6*time_axis,netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/CT_DEER_mix_28_36/CT_DEER_mix_28_36');
deernet(deer_trace,1e-6*time_axis,netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/deer_252cl_113scans/deer_252cl_113scans');
deernet(deer_trace,1e-6*time_axis,netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/deer_bi_19_50K/deer_bi_19_50K');
deernet(deer_trace,1e-6*time_axis,netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/deer_bi_28_50K/deer_bi_28_50K');
deernet(deer_trace,1e-6*time_axis,netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/deer_bi_36_50K/deer_bi_36_50K');
deernet(deer_trace,1e-6*time_axis,netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/deer_bi_50_50K/deer_bi_50_50K');
deernet(deer_trace,1e-6*time_axis,netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/deer_bi_oligo_n8_50K/deer_bi_oligo_n8_50K');
deernet(deer_trace(1:end/2),1e-6*time_axis(1:end/2),netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/deer_bi_oligo_n10_50K/deer_bi_oligo_n10_50K');
deernet(deer_trace(1:end/2),1e-6*time_axis(1:end/2),netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/deer_mixture_80K_m_8_5scans/deer_mixture_80K_m_8_5scans');
deernet(deer_trace,1e-6*time_axis,netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/deer_tri_36_50K/deer_tri_36_50K');
deernet(deer_trace,1e-6*time_axis,netset_location); drawnow();

[deer_trace,time_axis]=elexsys2deernet('sample_set_deeran/dOTP_5nm/dOTP_5nm');
deernet(deer_trace,1e-6*time_axis,netset_location); drawnow();

end

