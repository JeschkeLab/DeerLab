function [pass,maxerr] = test(opt)

% Check error control of regoperator() towards wrong inputs

N = 20;

% Pass 1: first input passed as a matrix
try
  regoperator(rand(100,10),2);
  pass(1) = false;
catch 
  pass(1) = true;
end

% Pass 2: order exceeds 3
try
  regoperator(N,4);
  pass(2) = false;
catch 
  pass(2) = true;
end

% Pass 3: no order specified
try
  regoperator(N);
  pass(3) = false;
catch 
  pass(3) = true;
end

pass = all(pass);

maxerr = NaN;
 

end