function window = kaiser_(N,beta)

narginchk(1,2);
if nargin < 2 || isempty(beta) 
    beta = 0.500;
end
N = round(N);
bessel = abs(besseli(0,beta));
odd = rem(N,2);
xind = (N-1)^2;
pos = fix((N+1)/2);
axis = (0:pos-1) + .5*(1-odd);
axis = 4*axis.^2;
window = besseli(0,beta*sqrt(1-axis/xind))/bessel;
window = abs([window(pos: - 1:odd+1) window])';

end