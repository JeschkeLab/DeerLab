
function [x,resnorm2,f,exitflag,Output,cov,sigx]=nlsqbnd(func,x,bl,bu,options,varargin)
%nlsqbnd solves non-linear least squares problems.
%   nlsqbnd attempts to solve problems of the form:
%   min  sum {func(x).^2}    where x and the values returned by func are vectors  
%    x                      
%
%   nlsqbnd has the same input parameters as lsqnonlin
%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   x = nlsqbnd(func,x0) starts at x0 and finds a minimum x to 
%   the sum of squares of the functions in func. func accepts input x 
%   and returns a vector of function values f evaluated at x.
%   NOTE: func should return func(x) and not the sum-of-squares 
%   sum(func(x).^2)). (func(x) is summed and squared implicitly in the
%   algorithm.) 
%
%   x = nlsqbnd (func,x0,lb,ub) defines a set of lower and upper bounds on
%   the design variables, x, so that the solution is in the range lb <= x<= ub.
%   Use empty vectors for lb and ub if no bounds exist. 
%   Set lb(i)= -Inf if x(i) is unbounded below; 
%   Set ub(i) = Inf if x(i) is unbounded above.
%
%   x = nlsqbnd (func,x0,lb,ub,options)  minimizes with the default
%   optimization parameters replaced by values in the structure options, an
%   argument created with the Matlab optimset function. 
%   Used options are (default values {  } :
%   Jacobian   =  [ on | {off} ]            (matlab standard)
%   'TolX'         =   {1e-6}                       
%   'TolRel'      =  {1e-6}                     (specific nlsqbnd)
%   'TolAbs'    =   {1e-6}                     (specific nlsqbnd)
%   'TolRank ',=  {1e-8}                      (specific nlsqbnd)
%    'MaxIter'   =   20*n                         (matlab standard)
%   'scaling'    =   [ {on} | off ]             (specific nlsqbnd)
%    'Newton'  =   [ {on} | off ]             (specific nlsqbnd)
%    'Display'  =    [{'off'}|'iter'|'final']  (matlab standard) 
%                             ('notify' not supported and starting point is
%                             also printed with the 'final' option.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%     Construction of func
%     **************************
%     func can be specified using @:
%        x = nlsqbnd (@myfun,x0)
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x) 
%             or [F,J] = myfun(x) when the 
%                              analytical jacobian J  is provided
%       F=.........
%       J=..........
%
%   FUN can also be an anonymous function:
%
%       x = lsqnonlin(@(x) sin(3*x),x0)
%
%   If func is parameterized, you can use anonymous functions to capture the 
%   problem-dependent parameters. Suppose you want to solve the non-linear 
%   least squares problem given in the function myfun, which is 
%   parameterized by its second argument c. Here myfun is an M-file 
%   function such as
%
%       function F = myfun(x,c)
%       F = [ 2*x(1) - exp(c*x(1))
%             -x(1) - exp(c*x(2))
%             x(1) - x(2) ];
%
%   To solve the least squares problem for a specific value of c, first 
%   assign the value to c. Then create a one-argument anonymous function 
%   that captures that value of c and calls myfun with two arguments. 
%   Finally, pass this anonymous function to nlsqbnd :
%
%       c = -1; % define parameter first
%       x = nlsqbnd (@(x) myfun(x,c),[1;1])
%
%    You can also use additionnal parameters with nlsqbnd 
%    x = nlsqbnd (@myfun,x0,lb,ub,options,P1, P2,...)
%    with a function such as F = myfun(x,P1,P2,...)
%
%    REMARK.- If the Jacobian J is provided analytically and can be computed
%    independently of f, it is efficient to compute J only if func is
%    called with two arguments.
%
%     Other nlsqbnd output arguments
%     [x,resnorm2,f,exitflag,output,cov,sigx] = nlsqbnd(func,x0,...) 
%     resnorm2 = ||f(x)||^2
%     f=f(x) at x;
%     exitflag describes the exit condition as follows :
%           the convergence criteria are :
%               1) relative predicted reduction in the objective function 
%                    is less than TolRel**2 
%               2) the sum of squares is less than TolAbs**2 
%               3) the relative change in x is less than TolX 
%               4) we are computing at noise level 
%                  the last digit in the convergence code (see below) indicates 
%                   how the last steps were computed 
%                   = 0 no trouble (gauss-newton the last 3 steps) 
%                   = 1 prank<>n at the termination point 
%                   = 2 the method of newton was used (at least) in the last step 
%                   = 3 the 2:nd but last step was subspace minimization but the 
%                       last two were gauss-newton steps 
%                   = 4 the steplength was not unit in both the last two steps 
%                the abnormal termination criteria are 
%                5) no. of iterations has exceeded maximum allowed iterations 
%                6) the hessian emanating from 2:nd order method is not pos def 
%                7) the algorithm would like to use 2:nd derivatives but is 
%                    not allowed to do that 
%                8) an undamped step with newtons method is a failure 
%                9) the latest search direction computed using subspace 
%                    minimization was not a descent direction (probably caused 
%                    by wrongly computed jacobian) 
%                exitflag integer scalar that indicate why the return is taken 
%                =10000  convergence due to criterion no. 1 
%                = 2000  convergence due to criterion no. 2 
%                =  300  convergence due to criterion no. 3 
%                =   40  convergence due to criterion no. 4 
%                =    x   where x equals 0,1,2,3 or 4 
%
%               <0   indicates that no convergence criterion is fulfilled 
%                       but some abnormal termination criterion is satisfied 
%               =  -1  if m<n or n<=0 or m<=0 or mdc<m or mdg<n or max<=0 
%                     or scale<0 or tol<0 or any of epsilon values <0 
%                     or invalid starting point   on entry 
%                =  -2   termination due to criterion no. 5 
%                =  -3   termination due to criterion no. 6 
%                =  -4   termination due to criterion no. 7 
%                =  -5   termination due to criterion no. 8 
%                =  -6   termination due to criterion no. 9 
%                =  -7   there is only one feasible point, namely 
%                           x(i)=bl(i)=bu(i)  i=1,2,....,n 
%                =  -20  termination due to user stop indicator f(x) not
%                           computable
%                =  -30  termination due to user stop indicator J(x) not
%                           computable
%
%      ouput is structure with
%           Output.iterations  : number of iterations
%           Output.funcCount : total number of function evaluations 
%           Output.Jrank : rank of J at optimum
%           Output.fCountJac : number of function evaluations 
%                   caused by Jacobian estimation by difference method 
%           Output.fCountHes : number of function evaluations
%                   caused by using 2:nd derivative information
%           Output.fLinesrh : number of function evaluations
%                    caused by the linesearch algorithm
%           Output.evalJ : number of analytical Jacobian evaluations
%      cov  contains inv(J'*J) at optimum
%      sigx is the standard deviation of x computed as
%                   sqrt(diag(resnorm*cov/(length(f)-length(x))))
%
%       Author : Alain Barraud 2009
%
%       This code is based upon a slightly modified version of ELSUNC
%        written by  Lindstrom and Wedin
%        Institute of information processing university of UME, SWEDEN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(x);
%options check
defaultopt = struct('Jacobian','off','TolX',1e-6,'TolRel',1e-6,'TolAbs',1e-6,'TolRank',...
    1e-8,'MaxIter',20*n,'scaling','on','Newton','on','Display','off');
if (~exist('options','var')) 
    options=defaultopt;
else
    f = fieldnames(defaultopt);
    for i=1:length(f),
        if (~isfield(options,f{i})||(isempty(options.(f{i})))), 
            options.(f{i})=defaultopt.(f{i}); end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=x(:);
if nargin<4,bu=repmat(inf,n,1);end
if isempty(bu),bu=repmat(inf,n,1);end
if nargin<3,bl=-repmat(inf,n,1);end
if isempty(bl),bl=-repmat(inf,n,1);end
bl=bl(:);bu=bu(:);
I=find(bl>=bu);if ~isempty(I),error('bl must be < bu');end
I=find(x>bu);if ~isempty(I),x(I)=bu(I);end
I=find(x<bl);if ~isempty(I),x(I)=bl(I);end
F=func(x,varargin{:});m=length(F);
if m==0,error('function is not computable at x0');end
%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(options.Display,'off')
    options.iprint=0;
elseif strcmp(options.Display,'iter')
    options.iprint=1;
elseif strcmp(options.Display,'final')  
    options.iprint=-1;
end
param=[options.TolRank;options.TolRel;options.TolAbs;options.TolX];
param(5:7)=[options.iprint;1;options.MaxIter];
if strcmp(options.Jacobian,'on')
    Joption=1;
else
    Joption=0;
end
if strcmp(options.Newton,'on')
    param(8)=-1;
else
    param(8)=1;
end
if strcmp(options.scaling,'on')
    param(9)=1;
else
    param(9)=0;
end
if param(8)<0%newton step available
    mdw = n * n + 5 * n + 3 * m + 6;
else
    mdw = 6 * n + 3 * m + 6;
end
mdp= 11 + 2 * n;
nc=max(n,4);
evalJ=0;
W=zeros(mdw,1);P=zeros(mdp,1);C=zeros(m,nc);
param(10)=evalJ;param(11)=Joption;
W(1:11)=param;
[x,f,exitflag,output,cov]=elsunc(func,x,bl,bu,W,C,varargin{:});
resnorm2=2*output(7);
if m==n
    sigx=repmat(inf,n,1);
else
    sigx=sqrt(diag(cov)*resnorm2/(m-n));
end
Output.iterations=output(1);
Output.funcCount=output(2);
Output.Jrank=output(6);
Output.fCountJac=output(3);
Output.fCountHes=output(4);
Output.fLinesrh=output(5);
Output.evalJ=output(9);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
