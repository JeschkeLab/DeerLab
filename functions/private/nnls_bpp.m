% Non-Negative Least Squares using Block Principal Pivoting
%
% Input:
%        AtA    nxn matrix
%        AtB    nx1 vector or nxk matrix
%        x0     initial x, nx1 vector or nxk matrix (default all-zero)
%
% Output:
%        x      solution vector(s), nx1 or nxk

% References:
% - Portugal, Judice, Vicente, Mathematics of Computation, 1994, 63, 625-643
%   A comparison of block pivoting and interior-point algorithms for linear
%   least squares problems with nonnegative variables
%   https://doi.org/10.1090/S0025-5718-1994-1250776-4
% - Kim, Park,SIAM J. Sci. Comput. 2011, 33(6), 3261-3281
%   Fast Nonnegative Matrix Factorization: An Active-Set-Like Method and Comparisons
%   https://doi.org/10.1137/110821172

function x = nnls_bpp(AtA,AtB,x0)

%Turn off warnings to avoid ill-conditioned warnings 
warning('off','all')

% Size checks
%-------------------------------------------------------------------------------
[n1,n2] = size(AtA);
if n1~=n2
  error('AtA must be a square matrix. You gave a %dx%d matrix.',n1,n2);
end
[n,k] = size(AtB);
if n~=n1
  error('AtB must have the same number of rows as AtA. You gave %d instead of %d',n,n1);
end
if nargin>2
  [n5,n6] = size(x0);
  if (n5~=n) || (n6~=k)
    error('x0 must have the same size as AtB (%dx%d). You gave a %dx%d array.',...
      n,k,n5,n6);
  end
end

% Loop over multiple right-hand sides
%-------------------------------------------------------------------------------
if k > 1
  for k_ = k:-1:1
    x(:,k_) = nnls_bpp(AtA,AtB(:,k_),x0(:,k_));
  end
  return
end

% Calculate initial solution
%-------------------------------------------------------------------------------
x = zeros(n,1);
if nargin < 3
  Fset = false(n,1);
  y = -AtB;
else
  Fset = x0 > 0;
  x = zeros(n,1);
  x(Fset) = AtA(Fset,Fset)\AtB(Fset);
  y = AtA*x - AtB;
end

% Determine infeasible variables ( = variables with negative values)
xFnegative = (x < 0) & Fset;
yGnegative = (y < 0) & ~Fset;
nInfeasible = sum(xFnegative | yGnegative);

% Iterative algorithm
%-------------------------------------------------------------------------------
p = 3;
t = inf;
while nInfeasible > 0 % iterate until no infeasible variables left
  
  % Swap full blocks in/out of F set; or swap a single element as backup plan.
  if nInfeasible < t
    t = nInfeasible;
    p = 3;
    Fset(xFnegative) = false;
    Fset(yGnegative) = true;
  else
    if p >= 1
      p = p - 1;
      Fset(xFnegative) = false;
      Fset(yGnegative) = true;
    else
      idx_ = find(xFnegative | yGnegative,1,'last');
      Fset(idx_) = ~Fset(idx_);
    end
  end
  
  % Solve linear system over F set
  x = zeros(n,1);
  x(Fset) = AtA(Fset,Fset)\AtB(Fset);
  y = AtA*x - AtB;
  
  % Determine infeasible variables
  xFnegative = (x < 0) & Fset;
  yGnegative = (y < 0) & ~Fset;
  nInfeasible = sum(xFnegative | yGnegative);
  
  
  %Turn off warnings to avoid ill-conditioned warnings
  warning('on','all')

end
