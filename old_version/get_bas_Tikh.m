function [K,t,r] = get_bas_Tikh(N,tmin,tmax,rmin,rmax)
% [K,t,r] = get_bas_Tikh(dimension,dt,rmin,rmax)
%
% Compute dipolar kernel matrix for use with Tikhonov regularization.
%
% Input:
%   N        kernel dimension (N x N)
%   tmin     minimum time in 탎, defaults to 0 microseconds
%   tmax     maximum time in 탎, defaults to (N-1)*0.008 microseconds
%   rmin     minimum distance in nm, defaults to 1.5 nm
%   rmax     maximum distance in nm, defaults to 10 nm
%
% Output:
%   K        kernel matrix, size N x N
%   t        time vector, microseconds
%   r        distance vector, nm

persistent cache
cacheSize = 10;

if nargin<1
  error('Kernel matrix size required.');
end

if ~exist('tmin','var') || isempty(tmin)
  tmin = 0;
end

if ~exist('tmax','var') || isempty(tmax)
  tmax = (N-1)*0.008; % 탎
end

dt = (tmax-tmin)/(N-1); % time increment, in 탎
nu_dip = 52.0410; % dipolar frequency at 1 nm for g=ge, MHz

if ~exist('rmin','var') || isempty(rmin)
  % minimum detectable distance, the largest dipolar frequency (shoulders of
  % the Pake doublet is 85% of the Nyquist frequency)
  rmin = (4*dt*nu_dip/0.85)^(1/3);
end

if ~exist('rmax','var') || isempty(rmax)
  % maximum detectable distance, 2 microseconds dipolar evolution time
  % corresponds to 6 nm
  rmax = 6*(N*dt/2)^(1/3);
end

if ~isempty(cache)
  for k = 1:numel(cache)
    c = cache(k);
    if all([N tmin tmax rmin rmax]==c.parameters)
      K = c.K;
      t = c.t;
      r = c.r;
      return
    end
  end
end

nt = N; % length of time axis
nr = N; % length of distance axis
t = linspace(tmin,tmax,nt).'; % time axis, in 탎
r = linspace(rmin,rmax,nr).'; % distance axis in nm
K = zeros(nt,nr); % kernel matrix
wdd = 2*pi*nu_dip./r.^3; % perpendicular dipolar angular frequencyies for all r

kernelCalcMethod = 2;
switch kernelCalcMethod
  case 1
    % Method using explicit numerical powder average (slow)
    %----------------------------------------------------------
    ntheta = 1001;
    costheta = linspace(0,1,ntheta);
    for ir = 1:nr
      Kr = 0;
      for l = 1:ntheta % average over cos(theta) angle (powder average)
        ww = wdd(ir)*(1-3*costheta(l)^2); % dipolar frequency at current theta
        Kr = Kr + cos(ww*t); % add kernel contribution
      end
      % normalize dipolar time evolution trace
      K(:,ir) = Kr/Kr(1);
    end
  case 2
    % Method using Fresnel integrals (fast)
    %----------------------------------------------------------
    % using John D'Errico's implementation of the Fresnel integrals
    % from the Matlab file exchange
    % Equations: see Edwards/Stoll, J.Magn.Reson. 270, 87-97 (2016), Eq.(6)
    % https://doi.org/10.1016/j.jmr.2016.06.021
    for ir = 1:nr
      ph = wdd(ir)*abs(t);
      y = sqrt(6*ph/pi);
      K(:,ir) = (cos(ph).*fresnelC(y)+sin(ph).*fresnelS(y))./y; % div by zero for y=0
    end
    % Correct at t = 0, since the Fresnel expression
    % involves division by zero for t = 0.
    K(t==0,:) = 1;
  otherwise
    error('Unknown calculation method for dipolar kernel.')
end

% Update cache
%-------------------------------------------------------------------------------
% If cache is full, drop oldest cached result.
if numel(cache)>=cacheSize
  cache(1) = [];
end

% Store result in cache.
idx = numel(cache)+1;
cache(idx).parameters = [N tmin tmax rmin rmax];
cache(idx).K = K;
cache(idx).t = t;
cache(idx).r = r;
