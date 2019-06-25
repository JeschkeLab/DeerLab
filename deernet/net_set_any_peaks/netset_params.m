% DEER library parameters
parameters.min_dist=15;      % Minimum distance, Angstrom
parameters.max_dist=50;      % Maximum distance, Angstrom
parameters.max_time=2e-6;    % DEER trace length
parameters.ndistmax=3;       % Maximum number of distances
parameters.npoints=256;      % Data vector size
parameters.noise_lvl=0.05;   % RMS noise, fraction of MD
parameters.min_exch=0.0;     % Maximum exchange coupling
parameters.max_exch=0.0;     % Minimum exchange coupling
parameters.min_bdim=2.0;     % Minimum background dimensionality
parameters.max_bdim=3.5;     % Maximum background dimensionality
parameters.min_fwhm=0.05;    % Minimum peak width, fraction of distance
parameters.max_fwhm=0.50;    % Maximum peak width, fraction of distance
parameters.min_skew=-3.0;    % Minimum peak skew, fraction of distance
parameters.max_skew=+3.0;    % Maximum peak skew, fraction of distance
parameters.max_mdep=0.6;     % Maximum modulation depth
parameters.min_mdep=0.1;     % Minimum modulation depth
parameters.min_brate=0;      % Minimum bg decay rate, s^-1
parameters.max_brate=0.5e6;  % Maximum bg decay rate, s^-1

% Network parameters
parameters.layer_sizes=[256 256 256 256];
parameters.lastlayer='logsig';
parameters.method='trainscg';

