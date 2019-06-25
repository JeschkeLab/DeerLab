function [nu,spc,window] = dft(t,signal,options)
% function [nu,spc,window] = dft(t,signal,options)
%
% Mathematically correct 1D Fourier transformation of a time-domain signal,
% providing frequency axis and spectrum
%
% t         time axis, if argument is empty, axis by point numbers
% signal    input signal, can be real or complex, for real signals, the
%           spectrum is returned only at positive frequencies, this can be
%           changed by options
%           signal should have an even number of points
% options   defines pre-processing, optional, all in-existent fields are
%           assigned default values
%           .force_complete     if true, spectrum at positive and negative
%                               frequencies is returned, even if imaginary
%                               component of signal is negligible, defaults
%                               to false
%           .correct_offset     true, automatic offset correction, defaults
%                               to false
%           .zerofilling        if false, no zero filling is applied
%                               if true (1), zero filling to twice the 
%                               number of original points is performed, 
%                               if larger than number of signal points, 
%                               zero filling to this number, 
%                               defaults to true (1)
%           .apodization        apodization window, string
%                               'none' no apodization, default
%                               'hamming' hamming apodization
%                               'dolph-chebyshev' Dolph-Chebyshev window
%                               'kaiser' Kaiser window
%                               'lorentz-gauss' Lorentz-Gauss
%                               transformation
%                               'sinebell' (shifted) sinebell apodization
%           .ripple             sidelobe attenuation in dB for
%                               Dolph-Chebyshev window, defaults to 50 dB
%           .alpha              alpha parameter for Kaiser window
%           .decayed            time axis fraction, at which the
%                               apodization window ends, defaults to 1
%           .tau                time constant for Lorentz part of
%                               Lorentz-Gauss transformation, defaults to
%                               1/(3*t_max)
%           .sigma              width parameter for Gauss part of
%                               Lorentz-Gauss transformation, defaults to
%                               1/(sqrt(2 ln 2) tau)
%           .phase              in rad, phase shift for sinebell apodization,
%                               defaults to pi/6 (30 degree)
%           .mode               determines, what data are returned
%                               'absorption' absorption sepctrum, real part
%                               'real'       synonym to 'absorption'
%                               'dispersion' dispersion spectrum, imaginary
%                                            part
%                               'imaginary'  synonym to 'dispersion'
%                               'magnitude'  magnitude spectrum, absolute
%                                            value
%                               'absolute'   synonym to 'magnitude'
%                               'power'      power spectrum
%                               'complex'    complex spectrum
%                               defaults to 'complex'
%
% nu        frequency axis, unit is the inverse time unit, row vector
% spc       spectrum, complex row vector
% window    window function used in apodization, if any, otherwise ones
%
% G. Jeschke, 2011

tol=1e-6; % tolerance for negligible imaginary input components

default_ripple=50; % default ripple attenuation for Dolph-Chebyshev window
default_alpha=2; % default alpha parameter of Kaiser apodization window
default_phase=10*pi/180; % default phase shift for sinebell apodization

nu=[];
spc=[];
window=ones(size(signal));

% return empty output, if signal input is empty
if isempty(signal),
    return
end;

% return empty output, if signal is not a vector
[m,n]=size(signal);
if n==1 && m>1, % make row vector, if column vector is given
    signal=signal.';
elseif m>1,
    fprintf(2,'Warning: fft_true called with 2D signal. Returning empty output.');
    return
end;

% define time axis, if none is given
if isempty(t),
    t=0:1:length(signal)-1;
end;

% determine, if the signal has significant imaginary components
complex_ratio=sum(abs(imag(signal)))/sum(abs(real(signal)));

% set default options
if nargin<3,
    options.force_complete=false;
end;

if ~isfield(options,'force_complete'),
    options.force_complete=false;
end;
if ~isfield(options,'correct_offset'),
    options.correct_offset=false;
end;
if ~isfield(options,'zerofilling'),
    options.zerofilling=1;
end;
if ~isfield(options,'apodization'),
    options.apodization='none';
end;
if ~isfield(options,'decayed'),
    options.decayed=1;    
end;
if ~isfield(options,'mode'),
    options.mode='complex';    
end;

% determine, if complete spectrum should be returned
if options.force_complete || complex_ratio>tol,
    complete=true;
else
    complete=false;
end;

% determine, how many zeros should be added
zf=length(signal);
if options.zerofilling,
    if options.zerofilling>length(signal),
        zf=options.zerofilling;
    else
        zf=2*length(signal);
    end;
end;
ntmax=t(1)+(t(end)-t(1))*(zf-1)/(length(signal)-1);
tzf=linspace(t(1),ntmax,zf);
zf=zf-length(signal);


% scale first point by 1/2 (formal definition of DFT leads to wrong integral of the spectrum)
signal(1)=signal(1)/2;

% do the optional apodization
switch lower(options.apodization),
    case 'none'
    case 'hamming'
        [signal,window]=private_hamming(signal,options.decayed);
    case 'dolph-chebyshev'
        if ~isfield(options,'ripple'),
            options.ripple=default_ripple;
        end;
        [signal,window]=private_dolph_chebyshev(signal,options.ripple);
    case 'kaiser'
        if ~isfield(options,'alpha'),
            options.alpha=default_alpha;
        end;
        [signal,window]=private_kaiser(signal,options.alpha);
    case 'lorentz-gauss'
        if ~isfield(options,'tau'),
            options.tau=max(t)/3;
        end;
        if ~isfield(options,'sigma'),
            options.sigma=1/(sqrt(2*log(2))*options.tau);
        end;
        [signal,window]=private_Lorentz_Gauss(signal,t,options.tau,options.sigma);
    case 'sinebell'
        if ~isfield(options,'phase'),
            options.phase=default_phase;
        end;
        [signal,window]=private_sinebell(signal,t,options.phase);
    otherwise
        fprintf(2,'Warning: Unknown apodization window "%s" requested.\n',options.apodization);
end;

% do offset correction, if requested
if options.correct_offset,
    signal=signal-mean(signal);
end;

% do the zero filling
signal=[signal zeros(1,zf)];
t=tzf;

% and the fft
spc=fft(signal);

% determine frequency range and shift or truncate spectrum, if necessary
T=t(2)-t(1);
nu=0:1/(length(t)*T):1/T;
nu=nu(1:length(spc));
if complete,
    spc=fftshift(spc);
    nu=fftshift(nu);
    nu(nu>=(1/(2*T)))=nu(nu>=(1/(2*T)))-1/T;
else
    spc=spc(1:ceil(length(signal)/2));
    nu=nu(1:ceil(length(signal)/2));
end;


% convert to required output mode, if other than 'complex'

switch lower(options.mode)
    case {'absorption','real'}
        spc=real(spc);
    case {'dispersion','imaginary'}
        spc=imag(spc);
    case {'magnitude','absolute'}
        spc=abs(spc);
    case 'power'
        spc=abs(spc).^2;
end;

function [spcn,hamm] = private_hamming(spc,r)
%HAMMING	Hamming apodization of a spectrum, renamed private_hamming 
%           to avoid conflicts with the Signal Processing Toolbox 
%	spcn = private_hamming(spc,r)
%
%  r part at which hamming window has decayed to zero
%	

%	Gunnar Jeschke, 1998

[m,n]=size(spc);
spcn=spc;
rad=round(r*n);
arg=linspace(0,pi,rad);
hamm=0.54*ones(1,rad)+0.46*cos(arg);
if rad>n, hamm=hamm(1,1:n); end;
if rad<n, hamm=[hamm zeros(1,n-rad)]; end;
for k=1:m,
   spcn(k,:)=hamm.*spcn(k,:);
end;

function [spcn,window] = private_dolph_chebyshev(spc,r)
%	Dolph-Chebyshev apodization of a spectrum
%	[spcn,window] = private_dolph_chebyshev(spc,r)
%
%  r ripple-to-main lobe ratio in dB, defaults to 20
%	

%	Gunnar Jeschke, 2011

if ~exist('chebwin','file'),
    [~,n]=size(spc);
    window=ones(1,n);
    spcn=spc;
    fprintf(2,'Warning: Dolph-Chebyshev window not available.\n');
    fprintf(2,'No apodization performed.\n');
    fprintf(2,'(Probably signal processing toolbox is missing)\n');    
    return
end;

if nargin<2,
    r=default_ripple;
end;

[m,n]=size(spc);
spcn=spc;
window=fftshift(chebwin(2*n,r)');
window=window(1:n);
for k=1:m,
   spcn(k,:)=window.*spcn(k,:);
end;

function [spcn,window] = private_kaiser(spc,alpha)
%

[m,n]=size(spc);
spcn=spc;
arg=linspace(0,1,n);
arg=sqrt(1-arg.^2);
window=besseli(0,pi*alpha*arg)/besseli(0,pi*alpha);
for k=1:m,
   spcn(k,:)=window.*spcn(k,:);
end;

function [spcn,window] = private_Lorentz_Gauss(spc,t,tau,sigma)
%

if nargin<2,
    t=0:length(spc)-1;
end;
if nargin<3,
    tau=max(t)/3;
end;
if nargin<4,
    sigma=1/(sqrt(2*log(2))*tau);
end;
window=exp(t/tau-sigma^2*t.^2/2);

[m,~]=size(spc);
spcn=spc;
for k=1:m,
   spcn(k,:)=window.*spcn(k,:);
end;

function [spcn,window] = private_sinebell(spc,t,phase)
%

if nargin<2,
    t=0:length(spc)-1;
end;
if nargin<3,
    phase=default_phase;
end;
window=sin(pi*t/max(t)+phase);

[m,~]=size(spc);
spcn=spc;
for k=1:m,
   spcn(k,:)=window.*spcn(k,:);
end;