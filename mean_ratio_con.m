function [c,ceq] = mean_ratio_con(x,handles,r,lb,ub,exflag,mcon,kcon,ctol)
%Nonlinear constraints placed on mean distance ratio for any 2-distance 
%fit distibution
%
% Written by H.C. Hyde, 2010 

%Get lower and upper constraints from 'mcon' structure
lb=mcon.lb;
ub=mcon.ub;

%Get mean distance parameter indices (into 'x') from 'kcon' array
k1=kcon(1);  %index of <r1> 
k2=kcon(2);  %index of <r2>

if ~isempty(lb) & isempty(ub)
    casenum=1;
elseif isempty(lb) & ~isempty(ub) 
    casenum=2;
elseif ~isempty(lb) & ~isempty(ub)
    casenum=3;
else
    casenum=4;
end
    
%Set inequality constraints in matrix form: C(X)<=0
switch casenum
    case 1   %one constraint:   lower bound (lb)
        c = -(x(k2)/x(k1)-lb);    %constraint #1: r2/r1 >= lb
        
    case 2   %one constraint:   upper bound (ub)
        c =  (x(k2)/x(k1)-ub);    %constraint #2: r2/r1 <= ub 
        
    case 3   %two constraints: lower and upper bounds (lb,ub)
        c = [-(x(k2)/x(k1)-lb);   %constraint #1: r2/r1 >= lb
              (x(k2)/x(k1)-ub)];  %constraint #2: r2/r1 <= ub
          
    case 4  %no constraints
        c = [];
end

%Set equality constraints in matrix form: C(X)=0
ceq = [];