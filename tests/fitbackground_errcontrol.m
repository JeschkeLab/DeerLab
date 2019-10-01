function [err,data,maxerr] = test(opt,olddata)

t = linspace(0,10,100);
lam0 = 0.25;
k = 0.05;
V = dipolarsignal(t,3,'Background',td_exp(t,k),'moddepth',lam0);

tstart = 4.2424;
tend = 10.0000;

twrong = 5;
try
fitbackground(V,t,@td_exp,[tstart tend twrong]);
err(1) = true;
catch 
err(1) = false;    
end

model = rand(12,1);
try
fitbackground(V,t,model,tstart);
err(2) = true;
catch 
err(2) = false;    
end

try
fitbackground(V,t,@td_exp,'moddepth',[0.5 0.1]);
err(3) = true;
catch 
err(3) = false;    
end

try
fitbackground(V,t,@td_exp,'moddepth',2);
err(4) = true;
catch 
err(4) = false;    
end

try
fitbackground(V,t);
err(5) = true;
catch 
err(5) = false;    
end

try
fitbackground(V,t,@td_exp,[tend tstart]);
err(6) = true;
catch 
err(6) = false;    
end

err = any(err);
maxerr = 0;
data = [];

end