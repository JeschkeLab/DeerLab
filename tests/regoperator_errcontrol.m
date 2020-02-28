function [pass,maxerr] = test(opt)

N = 20;
try
  regoperator(rand(100,10),2);
  err(1) = true;
catch 
  err(1) = false;
end
try
  regoperator(N,4);
  err(2) = true;
catch 
  err(2) = false;
end
try
  regoperator(N);
  err(3) = true;
catch 
  err(3) = false;
end

pass = all(err);
maxerr = 0;
 

end