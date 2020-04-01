% Compute linear combination of 3D noncentral chi distrubtion
% (also known as 3D Rice distribution)

% 3D Rice distribution = non-central chi distribution

function P = multirice3d(r,nu,sig,a)

N = numel(nu);
if numel(sig)~=N
  error('Number of widths and number of centers must be equal.');
end
if numel(a)~=N
  error('Number of amplitudes and number of centers must be equal.');
end
if any(a<0)
  error('Amplitudes must not be negative.');
end

n = 3; % degrees of freedom

r = r(:);
P = 0;
for k = 1:N
    s2 = sig(k).^2;
    I_scaled = besseli(n/2-1,nu(k)*r/s2,1);
    P = P + a(k)*nu(k)^(n/2-1)./s2*r.^(n/2).*exp(-(r.^2+nu(k)^2)/(2*s2)+nu(k)*r/s2).*I_scaled;
end

P(P<0) = 0;
P = P(:);

% Normalize
if ~all(P==0)
  P = P/sum(P)/mean(diff(r));
end
