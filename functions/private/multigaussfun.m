% Calculate linear combination of Gaussian functions

function P = multigaussfun(r,r0,fwhm,a)

N = numel(r0);
if numel(fwhm)~=N
  error('Number of widths and number of centers must be equal.');
end
if numel(a)~=N
  error('Number of amplitudes and number of centers must be equal.');
end

r = r(:);
P = 0;
for k = 1:numel(r0)
    sig = fwhm(k)/sqrt(2*log(2));
    P = P + a(k)*sqrt(2/pi)*1/sig*exp(-2*((r-r0(k))/sig).^2);
end

% Normalize
if ~all(P==0)
  P = P/trapz(r,P);    
end
