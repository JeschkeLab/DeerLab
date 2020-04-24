function [pass,maxerr] = test(opt)

% Check that fitsignal() works with a dipolar evolution function

r = 2.5; % nm
excitewidth = 15; % MHz
t = linspace(0,2,501); % us

% (1) Reference signal with infinite bandwidth
Kinf = dipolarkernel(t,r,'ExcitationBandwidth',inf,'Method','grid');

% (2) Signal using Fresnel integrals and approx. bandwidth treatment
K1 = dipolarkernel(t,r,'ExcitationBandwidth',excitewidth,'Method','grid');
K2 = dipolarkernel(t,r,'ExcitationBandwidth',excitewidth,'Method','integral');

% (3) Manual calculation using numerical powder average (reference)
K0 = zeros(numel(t),numel(r));
nKnots = 10000;
costheta = linspace(0,1,nKnots);
wdd = 2*pi*52.04./r.^3; % Mrad/s
q = 1-3*costheta.^2;
for ir = 1:numel(wdd)
  D_ = 0;
  for itheta = 1:nKnots
    w = wdd(ir)*q(itheta);
    D_ = D_ + cos(w*abs(t)).*exp(-w.^2/excitewidth^2);
  end
  K0(:,ir) = D_/nKnots;
end
K0= K0/K0(1);

% Pass 1: limited bandwidth cases are not equal to infinit bandwidth case
pass(1) = ~isequal(Kinf,K0) & ~isequal(Kinf,K1) & ~isequal(Kinf,K2);
% Pass 2: computed kernels are equal to the reference
pass(2) = all(all(abs(K0-K1)<1e-4)) & all(all(abs(K0-K2)<1e-4));

pass = all(pass);

maxerr = NaN;

if opt.Display
    plot(t,K0,t,K1,t,K02,t,Kinf,'.k','LineWidth',1.5);
    legend('Fresnel','grid','integral','inf bw');
   xlabel('t [\mus]')
   ylabel('K(t,2.5nm)')
   grid on, axis tight, box on
   
end

end