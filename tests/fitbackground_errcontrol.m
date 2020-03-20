function [pass,maxerr] = test(opt)

% Check error control of fitbackground() towards wrong inputs

t = linspace(0,10,100);
lam0 = 0.25;
k = 0.05;
V = dipolarsignal(t,3,'Background',bg_exp(t,k),'moddepth',lam0);

tstart = 4.2424;
tend = 10.0000;

% Pass 1: one time too much
twrong = 5;
try
fitbackground(V,t,@bg_exp,[tstart tend twrong]);
pass(1) = false;
catch 
pass(1) = true;    
end

% Pass 2: not passing a valid model function handle
model = rand(12,1);
try
fitbackground(V,t,model,tstart);
pass(2) = false;
catch 
pass(2) = true;    
end

% Pass 3: passing several moodulation depths
try
fitbackground(V,t,@bg_exp,'moddepth',[0.5 0.1]);
pass(3) = false;
catch 
pass(3) = true;    
end

% Pass 4: passing a wrong modulation depth
try
fitbackground(V,t,@bg_exp,'moddepth',2);
pass(4) = false;
catch 
pass(4) = true;    
end

% Pass 5: not enough inputs
try
fitbackground(V,t);
pass(5) = false;
catch 
pass(5) = true;    
end

% Pass 6: non increasing vector
try
fitbackground(V,t,@bg_exp,[tend tstart]);
pass(6) = false;
catch 
pass(6) = true;    
end

pass = all(pass);
maxerr = NaN;
 

end