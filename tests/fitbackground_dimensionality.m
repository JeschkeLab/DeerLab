function [pass,maxerr] = test(opt)

% Check indifference of fitbackground() towards input dimensionality

t = linspace(0,5,100);
bckg = exp(-(0.5*t));

tstart = [t(1) t(end)];

[B1,~,param1,~] = fitbackground(bckg,t,@bg_strexp,tstart);
[B2,~,param2,~] = fitbackground(bckg.',t,@bg_strexp,tstart);
[B3,~,param3,~] = fitbackground(bckg,t.',@bg_strexp,tstart);
[B4,~,param4,~] = fitbackground(bckg.',t.',@bg_strexp,tstart.');

% Pass 1: all backgrounds are equal
pass(1) = isequal(B1,B2,B3,B4);
% Pass 1: all backgrounds are equal
pass(2) = iscolumn(B1) & iscolumn(B2) & iscolumn(B3) & iscolumn(B4);
% Pass 1: all parameter vectors are rows
pass(3) = ~iscolumn(param1) & ~iscolumn(param2) & ~iscolumn(param3) & ~iscolumn(param4);

pass = all(pass);

maxerr = NaN;
 

end