% Compute 3D noncentral chi distrubtion (aka 3D Rice distribution)

function P = rice3d(r,nu,sig)

n = 3; % degrees of freedom

I_scaled = besseli(n/2-1,nu*r/sig^2,1);
P = nu^(n/2-1)./(sig^2)*r.^(n/2).*exp(-(r.^2+nu^2)/(2*sig^2)+nu*r/sig^2).*I_scaled;

P(P<0) = 0;
P = P(:);
