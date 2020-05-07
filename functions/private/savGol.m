
function g = savGol (f, nl, nr, M)
% SAVGOL SavGol smoothes the data in the vector f by means of a
%        Savitzky-Golay smoothing filter.
%
%        g = savGol (f, nl, nr, M)
%        Input: f : noisy data
%        nl: number of points to the left of the reference point
%        nr: number of points to the right of the reference point
%        M : degree of the least squares polynomial
%
%        Output: g: smoothed data
%
%        Example:
%        g = savGol(f,16,16,4)
%
%        In many applications one is measuring a variable that is both 
%        slowly varying and corrupted by random noise. Then it is often 
%        desirable to apply a smoothing filter to the measured data in 
%        order to reconstruct the underlying smooth function. We assume 
%        that the noise is independent of the observed variable and that 
%        the noise obeys a normal distribution with zero mean and given 
%        variation.
%
%        W. H. Press and S. A. Teukolsky,
%        Savitzky-Golay Smoothing Filters,
%        Computers in Physics, 4 (1990), pp. 669-672.
%
if nargin==0
    disp('No input arguments.');
    help('savGol');
    return;
end
    
    
% matrix A
A = ones (nl+nr+1, M+1);
for j = M:-1:1,
  A(:, j) = (-nl:nr)' .* A (:, j+1);
end
% filter coefficients c
[Q, R] = qr (A,0);
c = Q(:, M+1) / R(M+1, M+1);
% smoothing of the noisy data
% Note that there are two equivalent ways to apply the Savitzky-Golay
% filter to the vector f.  In the first case we use a for-loop whereas
% in the second case we use the faster built-in function filter.
%
% g = f;
% n = size (f);
% for i = 1+nl:n-nr,
%   g (i) = c' * f (i-nl:i+nr);
% end
%
n = length(f);
g = filter(c (nl+nr+1:-1:1), 1, f);
g(1:nl) = f(1:nl);
g(nl+1:n-nr) = g(nl+nr+1:n);
g(n-nr+1:n) = f(n-nr+1:n);
