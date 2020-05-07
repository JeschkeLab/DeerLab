%  LEVENBERG-MARQUARD NON-LINEAR LEAST-SQUARES SOLVER
%
%   lmlsqnonlin is similar to lsqnonlin with the Levenberg-Marquardt algorithm
%
%   x = lmlsqnonlin(obj,x0);
%   x = lmlsqnonlin(obj,x0,lb);
%   x = lmlsqnonlin(obj,x0,lb,ub);
%   x = lmlsqnonlin(obj,x0,lb,ub,opt);
%   [x,resnorm,fval,exitflag,exargs,output,lambda,jacobian] = lmlsqnonlin(___)
%
%   where obj and x0 are the objective function and the inital guess.
%   lb and ub are the box constraints and can be left empty, or set to
%   infinity or -infinity, respectively, if no constraints are binding. 
%   opt allows to pass on other arguments as explained below.
%
%   Assembled from the code from the following projects:
%   
%   - Levenberg-Marquardt toolbox (version 3.0) by Alexander Dentler (BSD license)
%   - Jacobian toolbox (version 3.0) by Alexander Dentler (BSD license)
%
%   This function is a combination of both projects with function nesting,
%   code organization, and some refactoring.
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md). 
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [xCurrent,Resnorm,fval,exitflag,extra_arguments,output,lambda,Jx] = lmlsqnonlin(obj,x0,lb,ub,opt)

% Read input parameters
xForm = x0;
n = numel(x0);
x0 = x0(:);
t1 = cputime;

if exist('lb','var')
    if isempty(lb)
        lb = -inf(size(xForm));
    end
    lb(lb == -realmax) = -inf; 
else
    lb = -inf(size(xForm));
end

if exist('ub','var')
    if isempty(ub)
        ub = inf(size(xForm));
    end
        ub(ub == realmax) = inf; 
else
    ub = inf(size(xForm));
end

%===============================================
%  DEFAULT PARAMETERS
%===============================================

%--------------------------------------------------------------------------
%   Display parameters
%--------------------------------------------------------------------------
Display = 'none';             % different degrees of notifications, default gives final output
title = [];                   % user can add a title that will show up before the iteration (and renewals of the header) and on the figure name
IterDispRenewal = 30;         % renew header of iterative display after so many steps
OutputFcn = [];               % collects user specified output function which is either a function handle or a cell array with functions in each cell
PlotIterations = 0;           % create plots at each iteration
pltfn = [];                   % user can supply one function which takes x as an argument and plots it values
extra_arguments = cell(1,0);  % arguments created with each iteration that is passed as argument in the next iteration (for nested iterations)

%--------------------------------------------------------------------------
%   Convergence parameters
%--------------------------------------------------------------------------
% Iteration and function count
MaxIter = 100;               % maximal no of iterations
MaxFunEvals = 1000;           % maximal no of function calls
% Acceptance tolerances
AccTol = 0;                   % breakup optimization if Resnorm is below this value (accept function value)
FooTol = 1e-7;                % breakup optimization if gradients are below this value (first order optimality condition)
RelFooTol = 1e-6;             % breakup optimization if relative gradients are below this value (first order optimality condition)
IncrTol = 1e-6;               % new evaluated point shows enough improvement to be fully accepted and lambda does not decrease
TolFun = 1e-6;                 % breakup optimization if absolute Resnorm improvement falls below this level
RelTolFun = 1e-6;              % breakup optimization if relative Resnorm improvement falls below this level
TolX = 1e-6;                  % breakup optimization if all absolute changes in parameters fall below this level
RelTolX = 1e-6;               % breakup optimization if all relative changes in parameters fall below this level

%--------------------------------------------------------------------------
%   Levenberg-Marquardt parameters
%--------------------------------------------------------------------------
JacobianMethod = 'simple';    % user supplies jacobian as second output argument of fun if set to 'on', otherwise 'off' (default)
MaxStepNo = 1;                % no of steps taken to find jacobian when "Jacobiab" is set to 'limit' or 'extrapolation'
FinDiffRelStep = eps^(1/3);   % multiplies step size of non-supplied Jacobian finite difference step
TypicalX = 1;                 % scales step size of non-supplied Jacobian finite difference step
DerivativeCheck = 'off';      % if set to 'on' it checks Jacobian in first step
ScaleProblem = 'none';        % if set to 'Jacobian' the problem is rescaled, by 'none' its not.
InitDamping = 1e-2;           % initial dampening
MinDamping = 1e-7;            % minimal dampening
MaxDamping = 1e7;             % maximal dampening
FactDamping = 10;             % increases or decreases dampening in loop
MaxEigTol = 1e-6;             % if largest eigenvalue becomes smaller than this value we attempt to use contraction mapping
Broyden_updates = 'on';       % set to 'on' it gives Broyden updates for the Jacobian for every 2*n steps, set to 'off' it requires updates in each iteration,
conservative_updates = 1;     % set to 1 it will only enforce the tolerances for foo, stepsize or eigenvalue when we just updated Jacobian

%==================================================
%  DYNAMIC READ IN OF STRUCTURE THAT GIVES OPTIONS
%==================================================

fval=[];J=[];
if nargin>=5 % optional parameter structure opt has been provided
    if isstruct(opt)
        list=who;
        %  check specifically if bounds are also given in option structure
        %   and check if they are identical, otherwise give error
        if isfield(opt, 'lb')
            if ~isequal(opt.lb(:),lb(:)) && ~isequal(-inf(numel(x0),1),lb(:))
                error('Lower bounds given are not identical to lower bounds in option structure.')
            end
        end
        if isfield(opt, 'ub')
            if ~isequal(opt.ub(:),ub(:)) && ~isequal(inf(numel(x0),1),ub(:))
                error('Upper bounds given are not identical to upper bounds in option structure.')
            end
        end
        %  dynamic read in
        for i1=1:numel(list)
            if isfield(opt, char(list{i1}))
                eval(horzcat(matlab.lang.makeValidName(char(list(i1))),'=opt.',char(list{i1}),';'));
            end
        end
    end
end

%==================================================
%  ERROR CHECKING AND PREPARATION
%==================================================
% Boundary error check
if (numel(lb)>0 || numel(lb)>0) && (numel(x0)~=numel(lb) || numel(x0)~=numel(ub)) || any(lb(:)>=ub(:))
    error('Number of bounds does not correspond to number of variable to optimize, or bounds are reversed.')
end

lb=lb(:);
ub=ub(:);
if any(x0<lb) || any(x0>ub)
    warning('The guess is outside the specified domain. We correct for this. But don''t let that happen again.')
    x0(x0<lb)=lb(x0<lb)+eps;
    x0(x0>ub)=lb(x0>ub)-eps;
end
% Variable transformations
unbndguess = bnd2unbnd(x0);
transformback = @(y)unbnd2bnd(y);

% Objective
objunbdn=@(x,vara)obj(reshape(x,size(xForm)),vara{:});
objbnd=@(y,vara)obj(reshape(transformback(y),size(xForm)),vara{:});

% Jacobian matter
switch JacobianMethod
    case 'simple'
        Jacobian_method=1;
    case 'limit'
        Jacobian_method=2;
    case {'extrapolation','romberg'}
        Jacobian_method=3;
    case 'user-provided'
        Jacobian_method=4;
    otherwise
        error('Jacobian evaluation method not found.')
end
option_Jacobian.Jacobian_method=Jacobian_method;
option_Jacobian.FinDiffRelStep=FinDiffRelStep;
option_Jacobian.TypicalX=TypicalX;
option_Jacobian.MaxStepNo=MaxStepNo;
option_Jacobian.lb=lb;
option_Jacobian.ub=ub;
if ~isscalar(Broyden_updates)
    if ischar(Broyden_updates)
        switch Broyden_updates
            case 'off'
                Broyden_updates=0;
            otherwise
                Broyden_updates=2*n;
        end
    end
end
%   dampening
lambda=InitDamping;lambda_old=InitDamping;
%   display
switch Display
    case {'notify','notify-detailed'}
        prnt = 1;
    case {'none','off'}
        prnt = 0;
    case {'iter','iter-detailed'}
        prnt = 3;
    case {'final','final-detailed'}
        prnt = 2;
    case 'simplex'
        prnt = 4;
    otherwise
        prnt = 1;
end
option_Jacobian.dsp=2*(prnt==3)+(prnt==4);
prnt_fun=@(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,how)fprintf('%4.0f  %4.0f  % 8.2g   % 8.2g   % 8.2g   % 8.2g   % 8.2g      % 8.2g          % 8.2g     % 8.2g % 8.2g % 8.2g % 8.2g   %s  \n',x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,how);
txt=['\n                                                 Norm   Rel norm   First-order   Rel first-order      Largest                              Comment\n',...
    '  I#    F#      f(x)      Df(x)   relDf(x)    of step    of step    optimality        optimality   eigenvalue   lambda      rho    ratio          '];
header_wo_title=sprintf(txt);
if ~isempty(title)
    header = sprintf(horzcat('\n',title,txt));
else
    header = header_wo_title;
end
%==================================================
%  STARTING VALUES
%==================================================
flag_nochange = 1;
h=NaN(n,1);E=NaN(n,1);
iteration=0;funccount=1;Jacobian_counter=1;how='Initial evaluation';
xLMstep=NaN(n,1);yLMstep=NaN(n,1);newtonianstep=NaN(n,1);xgradient=NaN(n,1);
rho=NaN;ratio=NaN;fun_improv=NaN;rel_fun_improv=NaN;maxFoo=NaN;maxrelFoo=NaN;
maxAbsStep=NaN;maxAbsRelStep=NaN;MaxEigJJ=NaN;

%--------------------------------------------------------------------------
% First evaluation
%--------------------------------------------------------------------------
if ~isempty(fval) && (~isempty(J) || (isempty(J) && Jacobian_method<4))
    Resnorm=sum(fval.^2);
else
    [fval,Resnorm,J,extra_arguments]=eval_fun(unbndguess,extra_arguments);
end
option_Jacobian.m=numel(fval);

% Initial display
if prnt>1
    disp(header)
    prnt_fun(MaxIter,MaxFunEvals,AccTol,TolFun,RelTolFun,TolX,RelTolX,FooTol,RelFooTol,...
        MaxEigTol,MaxDamping,IncrTol,NaN,'Thresholds');
    prnt_first={iteration,funccount,Resnorm,fun_improv,rel_fun_improv,maxAbsStep,maxAbsRelStep,maxFoo,maxrelFoo,...
        MaxEigJJ,lambda,rho,ratio,how};
    prnt_fun(prnt_first{:});
end

%   good staring values?
if any(~isfinite(fval))
    disp('staring values: ')
    disp(num2str(unbnd2bnd(unbndguess)))
    disp('fval: ')
    disp(num2str(fval))
    error('Some starting value is not finite.')
end

%   derivative check
if strcmp(DerivativeCheck,'on') && Jacobian_method==4
    option_Jacobian_check=option_Jacobian;
    option_Jacobian_check.Jacobian_method=3;
    J_check=eval_Jacobian(@(x)objunbdn(x,extra_arguments),unbndguess,fval,option_Jacobian_check);
    dh=D_unbnd2bnd(unbndguess,lb,ub);
    err=bsxfun(@times,(J-J_check),1./dh);
    relerr=err./bsxfun(@plus,abs(transformback(unbndguess)),abs(transformback(unbndguess)'))/2;
    disp(horzcat('Derivative check gives largest error with ',...
        num2str(max(err(:))),' and largest relative error with ',num2str(max(relerr(:)))))
end

%==================================================
%  DYNAMIC LOOP TO OPTIMIZE
%==================================================

%   initial call of output functions and plotting function
OtptFcnVl.xInit=x0;OtptFcnVl.fvalInit=fval;
OtptFcnVl.pltfn=pltfn;OtptFcnVl.transformback=transformback;
OtptFcnVl.xForm=xForm;OtptFcnVl.title=title;
OtptFcnVl.ResnormHist=Resnorm;
OtptFcnVl.iteration=iteration;OtptFcnVl.funccount=funccount;
OtptFcnVl.Jacobian_counter=Jacobian_counter;OtptFcnVl.fval=fval;
OtptFcnVl.Resnorm=Resnorm;OtptFcnVl.J=J;OtptFcnVl.h=h;OtptFcnVl.E=E;
OtptFcnVl.xLMstep=xLMstep;OtptFcnVl.yLMstep=yLMstep;
OtptFcnVl.newtonianstep=newtonianstep;OtptFcnVl.xgradient=xgradient;
OtptFcnVl.extra_arguments=extra_arguments;
OtptFcnVl.rho=rho;OtptFcnVl.ratio=ratio;OtptFcnVl.lambda=lambda;
OtptFcnVl.fun_improv=fun_improv;OtptFcnVl.rel_fun_improv=rel_fun_improv;
OtptFcnVl.maxFoo=maxFoo;OtptFcnVl.maxrelFoo=maxrelFoo;
OtptFcnVl.maxabsstep=maxAbsStep;OtptFcnVl.maxabsrelstep=maxAbsRelStep;
OtptFcnVl.MaxEigJJ=MaxEigJJ;
% Callup functions
stop = eval_outputfun(OutputFcn,transformback(unbndguess),OtptFcnVl,'init');
if PlotIterations
    OtptFcnVl=eval_PlotIterations(transformback(unbndguess),OtptFcnVl,'init');
end
% Loop
howJ='user-supplied Jacobian';
while Resnorm>AccTol && funccount< MaxFunEvals && iteration < MaxIter && stop==false
    % Start new iteration
    iteration=iteration+1;
    maxAbsStep=NaN;maxAbsRelStep=NaN;MaxEigJJ=NaN;
    rho=NaN;ratio=NaN;fun_improv=NaN;rel_fun_improv=NaN;
    
    
    % Evaluate Jacobian at current point
    %---------------------------------------
    dh = D_unbnd2bnd(unbndguess);
    if isempty(J) && ~strcmp(JacobianMethod,'user-provided')
        [J,h,E,func_evals_Jacobian] = eval_Jacobian(@(x)objunbdn(x,extra_arguments),unbndguess,fval,option_Jacobian);
        funccount=funccount+func_evals_Jacobian;
        Jacobian_counter=1;
        howJ='full Jacobian update';
    end
    how=howJ;
    if prnt==4
        disp('Jacobian:')
        disp(num2str(bsxfun(@times,J,1./dh))) %     transform back to x space
        if Jacobian_method==2 || Jacobian_method==3
            disp('step size for Jacobian:')
            disp(num2str(h))
            disp('Error size:')
            disp(num2str(E))
        end
    end
    
    % Evaluate first order optimatlity
    %---------------------------------------
    xgradient = bsxfun(@times,-J,1./dh)'*fval;
    relFoo = xgradient./(1e-12+abs(transformback(unbndguess)));
    maxFoo = max(abs(xgradient));
    maxrelFoo = max(abs(relFoo));
    
    %If first order optimality fulfilled, break iteration
    if (maxFoo<FooTol || maxrelFoo<RelFooTol) && (Jacobian_counter==1 || ~conservative_updates)
        if maxFoo<FooTol && maxrelFoo<RelFooTol
            how='absolute and relative first order optimatility';
        elseif maxFoo<FooTol
            how='absolute first order optimatility';
        else
            how='relative first order optimatility';
        end
        break
    end
    
    % Set up the surrogate model
    JJ = J'*J;
    ygradient = -J'*fval;
    switch ScaleProblem
        case 'Jacobian'
            T=diag(diag(JJ));
        otherwise
            T=eye(size(JJ,1));
    end
    
    LMstepFcn = @(lambda)pinv((JJ+lambda*T),eps)*ygradient;
    MaxEigJJ = max(eig(JJ));
    if MaxEigJJ<MaxEigTol && (Jacobian_counter==1 || ~conservative_updates)
        how='largest eigenvalue';
        break
    end

    % Step size and stopping criteria for step size
    yLMstep = LMstepFcn(lambda);
    newtonianstep = LMstepFcn(0);
    xLMstep = transformback(unbndguess + yLMstep) - transformback(unbndguess);
    maxAbsStep = norm(xLMstep,2);
    maxAbsRelStep = norm(xLMstep./(1e-12+abs(transformback(unbndguess))),2);
    
    % If step size is too small, break iteration
    if (maxAbsStep<TolX || maxAbsRelStep<RelTolX) && (Jacobian_counter==1 || ~conservative_updates)
        if maxAbsStep<TolX && maxAbsRelStep<RelTolX
            how='absolute and relative step size';
        elseif maxAbsStep<TolX
            how='absolute step size';
        else
            how='relative step size';
        end
        break
    elseif (maxAbsStep<TolX || maxAbsRelStep<RelTolX)
        Jacobian_counter=Broyden_updates;
    end
    
    % Evaluate next point
    %---------------------------------------
    [LM_fval,LM_Resnorm,LM_J,LM_extra_arguments] = eval_fun(unbndguess+yLMstep,extra_arguments);
    funccount = funccount+1;
    if prnt==4
        disp('No, step size, New vs old parameters:')
        disp(num2str([ (1:n)' xLMstep transformback(unbndguess)+xLMstep transformback(unbndguess)]))
        disp(horzcat('Difference in norm of residuals :',num2str(LM_Resnorm-Resnorm)))
        disp(horzcat('New norm of residuals           :',num2str(LM_Resnorm)))
        disp(horzcat('Old norm of residuals           :',num2str(Resnorm)))
    end
    %  evaluate dampening
    fun_improv=LM_Resnorm-Resnorm;
    rel_fun_improv=LM_Resnorm/Resnorm-1;
    rho=(Resnorm-LM_Resnorm)/(2*yLMstep'*(lambda*yLMstep+ygradient));
    ratio=(Resnorm-sum((fval+J*yLMstep).^2))/(Resnorm-LM_Resnorm);
    lambda_old=lambda;
    if rho>IncrTol || TolFun>fun_improv || RelTolFun>rel_fun_improv
        %   note
        how = horzcat('*',how);
        %   good evaluation, decrease dampening
        lambda=max(lambda/FactDamping,MinDamping);
        flag_nochange=0;
        %  Jacobian update
        if Jacobian_method==4
            %  user supplied function also updated Jacobian
            J=LM_J;
            howJ='user-supplied Jacobian';
        else
            %  update jacobian or destroy it
            Jacobian_counter=Jacobian_counter+1;
            if Jacobian_counter>=Broyden_updates || isempty(J)%(2*n) || ~Broyden_updates
                J=[];
                howJ='full Jacobian update';
            else
                J=J+((LM_fval-fval-J*yLMstep)*yLMstep')/(yLMstep'*yLMstep);
                howJ='Broyden-type update';
            end
        end
        %  other values
        fval=LM_fval;
        Resnorm=LM_Resnorm;
        unbndguess=unbndguess+yLMstep;
        extra_arguments=LM_extra_arguments;
    elseif lambda==MaxDamping && (Jacobian_counter<=2 || ~conservative_updates)
        how='dampening';
        break
    else
        %   bad evaluation, increase dampening
        if Jacobian_counter>1 && Jacobian_method<4 && lambda==MaxDamping
            J=[];
            lambda=InitDamping;
            howJ='full Jacobian update';
        elseif Jacobian_counter>1 && Jacobian_method<4
            %   quick dampening as we use a Broyden updated jacobian which
            %   is quick, but suboptimal so we dont want to waste time with
            %   large steps that are imprecise.
            lambda=min(lambda*(FactDamping^2),MaxDamping);
            howJ='quick dampening';
        else
            lambda=min(lambda*FactDamping,MaxDamping);
            howJ='soft dampening';
        end
    end
    if Resnorm<=AccTol
        break
    end
    
    %  iterative display for LM
    if prnt>2
        if mod(iteration,IterDispRenewal)==0 || IterDispRenewal==1
            disp(header)
            prnt_fun(MaxIter,MaxFunEvals,AccTol,TolFun,RelTolFun,TolX,RelTolX,FooTol,RelFooTol,...
                MaxEigTol,MaxDamping,IncrTol,NaN,'Thresholds');
        end
        prnt_fun(iteration,funccount,Resnorm,fun_improv,rel_fun_improv,maxAbsStep,maxAbsRelStep,maxFoo,maxrelFoo,...
            MaxEigJJ,lambda_old,rho,ratio,how);
    end
    %  OUTPUT FUNCTIONS
    OtptFcnVl.ResnormHist=[OtptFcnVl.ResnormHist Resnorm];
    OtptFcnVl.iteration=iteration;OtptFcnVl.funccount=funccount;
    OtptFcnVl.Jacobian_counter=Jacobian_counter;OtptFcnVl.fval=fval;
    OtptFcnVl.Resnorm=Resnorm;OtptFcnVl.J=J;OtptFcnVl.h=h;OtptFcnVl.E=E;
    OtptFcnVl.xLMstep=xLMstep;OtptFcnVl.yLMstep=yLMstep;
    OtptFcnVl.newtonianstep=newtonianstep;OtptFcnVl.xgradient=xgradient;
    OtptFcnVl.rho=rho;OtptFcnVl.ratio=ratio;OtptFcnVl.lambda=lambda_old;
    OtptFcnVl.fun_improv=fun_improv;OtptFcnVl.rel_fun_improv=rel_fun_improv;
    OtptFcnVl.maxFoo=maxFoo;OtptFcnVl.maxrelFoo=maxrelFoo;
    OtptFcnVl.maxabsstep=maxAbsStep;OtptFcnVl.maxabsrelstep=maxAbsRelStep;
    OtptFcnVl.MaxEigJJ=MaxEigJJ;
    OtptFcnVl.extra_arguments=extra_arguments;
    stop=eval_outputfun(OutputFcn,transformback(unbndguess),OtptFcnVl,'iter');
    if PlotIterations
        OtptFcnVl=eval_PlotIterations(transformback(unbndguess),OtptFcnVl,'iter');
    end
    
end %End of main loop

%==================================================
%  FINAL ASSIGNMENT
%==================================================
xCurrent=reshape(transformback(unbndguess),size(xForm));
OtptFcnVl.iteration=iteration;OtptFcnVl.funccount=funccount;
if Resnorm<=AccTol
    how='FULL CONVERGENCE';
    exitflag=1;
elseif fun_improv>-TolFun || rel_fun_improv>-RelTolFun ...
        || maxFoo<FooTol || maxrelFoo<RelFooTol || maxAbsStep<TolX ...
        || maxAbsRelStep<RelTolX ||  MaxEigJJ<MaxEigTol ...
        || lambda==MaxDamping || rho<IncrTol
    exitflag=0;
elseif funccount==MaxFunEvals || iteration==MaxIter
    how='did not converge';
    exitflag=-1;
elseif stop
    how='Output function terminated algorithm';
    exitflag=1;
end
if nargout>5
    output=OtptFcnVl;
    output.how=how;
    output.nochange=flag_nochange;
    output.maxchange=max(abs(OtptFcnVl.xInit-xCurrent(:)));
    output.xCurrent=xCurrent;
    output.time_used_in_sec=cputime-t1;
    if Jacobian_method<4
        output.algorithm='levenberg-marquardt with Broyden rank-1 updates for Jacobian';
    else
        output.algorithm='levenberg-marquardt with user-supplied Jacobian';
    end
    output.ErrJacobian=E;
    output.StepsizeJacobian=h;
    output.extra_arguments=extra_arguments;
    if nargout>7
        if isempty(J) && Jacobian_method<4
            [J,h,E,func_evals_Jacobian]=eval_Jacobian(@(x)objunbdn(x,extra_arguments),unbndguess,fval,option_Jacobian);
            funccount=funccount+func_evals_Jacobian;
        end
        dh=D_unbnd2bnd(unbndguess,lb,ub);
        Jx=bsxfun(@times,J,1./dh);
        output.firstorder=(Jx'*fval);
        output.firstorderopt=max(abs((Jx'*fval)));
        output.Jacobian=Jx;
    end
end
%   close things
eval_outputfun(OutputFcn,transformback(unbndguess),OtptFcnVl,'done');
if PlotIterations
    eval_PlotIterations(transformback(unbndguess),OtptFcnVl,'done');
end
if prnt>0
    if prnt==3
        if iteration>1
            disp(header)
        end
        prnt_fun(MaxIter,MaxFunEvals,AccTol,TolFun,RelTolFun,TolX,RelTolX,FooTol,RelFooTol,...
            MaxEigTol,MaxDamping,IncrTol,NaN,'Thresholds');
        prnt_fun(prnt_first{:})
    elseif prnt==1
        disp(header_wo_title)
        %prnt_fun(prnt_first{:})
    end
    prnt_fun(iteration,funccount,Resnorm,fun_improv,rel_fun_improv,maxAbsStep,maxAbsRelStep,maxFoo,maxrelFoo,...
        MaxEigJJ,lambda_old,rho,ratio,how);
end


%==========================================================================
%                               SUBFUNCTIONS
%==========================================================================


    %==================================================
    % OBJECTIVE EVALUATION FUNCTION
    %==================================================
    function [fval,resnorm,J,args]=eval_fun(par0,args)
        J = [];
        % Evaluate the user-given objective function
        if strcmp(JacobianMethod,'user-provided')
            N = 1;
        else
            N = 0;
        end
        varargout = cell(1,1 + numel(args) + N);
        [varargout{:}] = objbnd(par0,args);
        fval = varargout{1};
        % Get Jacobian if provided by user
        if strcmp(JacobianMethod,'user-provided')
            J = varargout{2};
        end
        % If we got extra arguments...
        if numel(args)>0
            [args{:}] = deal(varargout{2+N:end});
        end
        
        % Get the l2-norm of the residuals
        resnorm=sum(fval.^2);
        
    end

    %==================================================
    % JACOBIAN EVALUATION FUNCTION
    %==================================================
    function [J,h,E,func_evals_Jacobian] = eval_Jacobian(ObjFcn,param,fval,opts)
        
        % Transform parameters to bounded space
        boundedGuess = unbnd2bnd(param);
        
        Jacobian_method = opts.Jacobian_method;
        if any(lb==boundedGuess) || any(ub==boundedGuess)
            Jacobian_method = 1;
        end
        opts.f_0 = fval;
        % Evaluate Jacobian
        switch JacobianMethod
            case 'simple'
                E=[];
                [J,h,func_evals_Jacobian]=jacobiansimple(ObjFcn,boundedGuess,opts);
            case 'limit'
                [J,h,func_evals_Jacobian,E]=jacobianlim(ObjFcn,boundedGuess,opts);
            case 'extrapolation'
                [J,h,func_evals_Jacobian,E]=jacobianext(ObjFcn,boundedGuess,opts);
        end
        
        % Transform back to y space
        df = D_unbnd2bnd(param);
        J = bsxfun(@times,J,df);
    end


    %==================================================
    % UNBOUND -> BOUND VARIABLE TRANSFORMATION
    %==================================================
    function X = unbnd2bnd(Y)
        % Transforms variable from -infinity to infinity to domain [lb,ub]
        % complements bnd2unbnd
        X=NaN(size(Y));
        for kk=1:numel(Y)
            if isfinite(lb(kk)) && isfinite(ub(kk)) % ...lower and upper bound
                X(kk) = (lb(kk)+ub(kk))/2+ (ub(kk)-lb(kk))/2*sin(2*Y(kk)/(ub(kk)-lb(kk)));
            elseif isfinite(lb(kk)) && ~isfinite(ub(kk)) % ...just lower bound
                X(kk)= lb(kk)-1 + sqrt(Y(kk).^2+1);
            elseif ~isfinite(lb(kk)) && isfinite(ub(kk)) % ...just upper bound
                X(kk)= ub(kk)+1 - sqrt(Y(kk).^2+1);
            else % ...no bounds
                X(kk)=Y(kk);
            end
        end
    end

    %==================================================
    % BOUND -> UNBOUND VARIABLE TRANSFORMATION
    %==================================================
    function    Y = bnd2unbnd(X)
        % Transforms variable from [lb,ub] to domain -infinity to infinity
        % complements unbnd2bnd
        Y=NaN(size(X));
        for jj=1:numel(X)
            if isfinite(lb(jj)) && isfinite(ub(jj)) % ...bounded on both ends
                Y(jj)=(ub(jj)-lb(jj))/2*asin((2*X(jj)-(ub(jj)+lb(jj)))/(ub(jj)-lb(jj)));
            elseif isfinite(lb(jj)) && ~isfinite(ub(jj)) % ...just lower bound
                Y(jj)=sqrt((lb(jj) -1 - X(jj)).^2-1);
            elseif ~isfinite(lb(jj)) && isfinite(ub(jj)) % ...just upper bound
                Y(jj)=sqrt((ub(jj) +1 - X(jj)).^2-1);
            else % ...no boundaries
                Y(jj)=X(jj);
            end
        end
    end



    function    df=D_unbnd2bnd(Y)
        df=NaN(1,numel(Y));
        for ii=1:numel(Y)
            if isfinite(lb(ii)) && isfinite(ub(ii)) % ...lower and upper bound
                df(ii) = cos(2*Y(ii)/(ub(ii)-lb(ii)));
            elseif isfinite(lb(ii)) && ~isfinite(ub(ii)) % ...just lower bound
                df(ii)= Y(ii)/sqrt(Y(ii).^2+1);
            elseif ~isfinite(lb(ii)) && isfinite(ub(ii)) % ...just upper bound
                df(ii)= -Y(ii)/sqrt(Y(ii).^2+1);
            else % ...no bounds
                df(ii)=1;
            end
        end
    end


end

%==========================================================================
%                           LOCAL FUNCTIONS
%==========================================================================

function    stop = eval_outputfun(OutputFcn,xGuess,optimValues,state)
stop=false;
if ~isempty(OutputFcn)
    if isa(OutputFcn,'function_handle')
        stop = OutputFcn(xGuess,optimValues,state);
    elseif isa(OutputFcn,'cell')
        stop=false(numel(OutputFcn),1);
        for i1=1:numel(OutputFcn)
            stop(i1)= OutputFcn{i1}(xGuess,optimValues,state);
        end
        stop=any(stop);
    end
end
end

function    OtptFcnVl = eval_PlotIterations(xGuess,OtptFcnVl,state)
xForm=OtptFcnVl.xForm;
switch state
    case 'init'
        fig=figure;
        %   optional figure name
        if ~isempty(OtptFcnVl.title)
            set(fig,'Name',OtptFcnVl.title)
        end
        subp(1)=subplot(3,2,1,'Parent',fig);title('current x');hold('on');
        bar(1:numel(xGuess),xGuess,'Parent',subp(1));
        
        subp(2)=subplot(3,2,2,'Parent',fig);title('current fval');hold('on');
        plot(1:numel(OtptFcnVl.fval),OtptFcnVl.fval,'Parent',subp(2));
        xlim(subp(2),[1 numel(OtptFcnVl.fval)])
        
        subp(3)=subplot(3,2,3,'Parent',fig);title('change from initial x');hold('on');
        subp(4)=subplot(3,2,4,'Parent',fig);title('change from initial fval');hold('on');
        subp(5)=subplot(3,2,5,'Parent',fig);title('attempted step');hold('on');
        
        if ~isempty(OtptFcnVl.pltfn)
            xForm(:)=xGuess;
            subp(6)=subplot(3,2,6,'Parent',fig);title('user-supplied');hold('on');
            temp=OtptFcnVl.pltfn(xForm);
            plot(temp,'Parent',subp(6));
            xlim(subp(6),[1 size(temp,1)])
        else
            subp(6)=subplot(3,2,4,'Parent',fig);title('history residual norm');hold('on');
        end
        OtptFcnVl.fig=fig;
        OtptFcnVl.subp=subp;
        pause(.001)
    case 'iter'
        for i1=1:numel(OtptFcnVl.subp)
            delete(get(OtptFcnVl.subp(i1), 'Children'));
        end
        %   parameter values
        bar(1:numel(xGuess),xGuess,'Parent',OtptFcnVl.subp(1));
        xlim(OtptFcnVl.subp(1),[.5 numel(xGuess)+.5])
        %   residuals
        plot(1:numel(OtptFcnVl.fval),OtptFcnVl.fval,'Parent',OtptFcnVl.subp(2));
        xlim(OtptFcnVl.subp(2),[1 numel(OtptFcnVl.fval)])
        %   change in parameter values from start
        bar(1:numel(xGuess),xGuess-OtptFcnVl.xInit,'Parent',OtptFcnVl.subp(3));
        xlim(OtptFcnVl.subp(3),[.5 numel(xGuess)+.5])
        %   change in residuals from start
        plot(1:numel(OtptFcnVl.fval),OtptFcnVl.fval-OtptFcnVl.fvalInit,'Parent',OtptFcnVl.subp(4));
        xlim(OtptFcnVl.subp(4),[1 numel(OtptFcnVl.fval)])
        %   attempted step
        bar(1:numel(xGuess),OtptFcnVl.xLMstep(:),'Parent',OtptFcnVl.subp(5));
        xlim(OtptFcnVl.subp(5),[.5 numel(xGuess)+.5])
        %   user supplied function
        if ~isempty(OtptFcnVl.pltfn)
            xForm(:)=xGuess;
            temp=OtptFcnVl.pltfn(xForm);
            plot(temp,'Parent',OtptFcnVl.subp(6));
            xlim(OtptFcnVl.subp(6),[1 size(temp,1)])
        else
            plot(OtptFcnVl.ResnormHist,'Parent',OtptFcnVl.subp(6));
        end
        OtptFcnVl.fig=OtptFcnVl.fig;
        OtptFcnVl.subp=OtptFcnVl.subp;
        pause(.001)
    case 'done'
        % Cleanup of plots, guis, or final plot
        close(OtptFcnVl.fig);
    otherwise
end
end


%==================================================
% JACOBIAN SIMPLE 
%==================================================
function [J,h,func_evals,f_0]=jacobiansimple(func,x,varargin)
%  DEFAULT PARAMETERS
f_0=[];                     %   value at evaluation point (can be user-supplied
lb=-inf(size(x));           %   lower boundary of domain of X
ub=inf(size(x));            %   upper boundary of domain of X
FinDiffRelStep=eps^(1/3);   %   relative differentiation step size
TypicalX=1;                 %   sets lower limit of differentiation step 
%                               size, similar to the optimization routines 
%                               in matlab 
%  DYNAMIC READ IN
if nargin>2    
    list=who;
    opt=[];arg_list=[];
    pos=1;flag_thereismore=1;
    %  check if argument gives us a structure
    if flag_thereismore
        if isstruct(varargin{pos})
            opt=varargin{pos};
            pos=pos+1;
            flag_thereismore=nargin>(pos+1);
        end
    end
    %  check for name-variable pairs
    if flag_thereismore
        if ((nargin-pos-1)/2)==fix((nargin-pos-1)/2)
            arg_list=varargin(pos:end);
        else
            error('No of arguments is off.')
        end
    end  
    %  add option structure variables if they are part of the list
    if ~isempty(opt)
        for i1=1:numel(list)
            if isfield(opt, char(list{i1}))
                eval(horzcat(matlab.lang.makeValidName(char(list(i1))),'=opt.',char(list{i1}),';'));
            end
        end
    end
    %  add name-value pair arguments if they are part of the list
    if ~isempty(arg_list)
        for i1=1:numel(arg_list)/2
            if ismember(arg_list{(i1-1)*2+1},list) && isnumeric(arg_list{(i1-1)*2+2})
                eval(horzcat(arg_list{(i1-1)*2+1},'=',num2str(arg_list{(i1-1)*2+2}),';'));
            end
        end
    end 
end
%  MAIN ROUTINE
% function value
func_evals=0;
if isempty(f_0)
    f_0=func(x);                   
    func_evals=func_evals+1;
end
%   size of independent
n=numel(x);  
%   size of dependent
m=numel(f_0);
%   size of increment before boundaries are considered
h=FinDiffRelStep.*max(abs(x),TypicalX);
%   allocate memory for the Jacobian matrix
J=zeros(m,n);                   
%  loop for each independent variable 
for k=1:n                       
    %  reference point calculation
    x1=x;                       
    %   boundary integrity
    if (ub(k)-x1(k))<h(k)         % feasibility
        if (x1(k)-lb(k))<h(k) 
            h(k)=-h(k);
        else
            if (ub(k)-x1(k))>=(x1(k)-lb(k))
                h(k)=ub(k)-x1(k);
            else
                h(k)=lb(k)-x1(k);
            end
        end
    end
    %   final increment in kth independent variable
    x1(k)=x1(k)+h(k);   
    %  step differentiation 
    J(:,k)=(func(x1)-f_0)/h(k);     
    func_evals=func_evals+1;
    h(k)=abs(h(k));
end
end

%==================================================
% JACOBIAN EXTRAPOLATED 
%==================================================
function [J,h,func_evals,E,S,D] = jacobianext(func,x,varargin)

%  SET OPTIONS
lb=-inf(size(x));           %   lower boundary of domain of X
ub=inf(size(x));            %   upper boundary of domain of X
FinDiffRelStep=eps^(1/3);   %   relative differentiation step size
TypicalX=1;                 %   sets lower limit of differentiation step 
%                               size, similar to the optimization routines 
%                               in matlab 
MaxStepNo=3;                %   max no of steps for the limiting process
toler=1e-6;                 %   absolute tolerances
%  DYNAMIC READ IN
D=[];   %   old Romberg matrix from earlier evaluation that can be reused. Inofficial part.
if nargin>2    
    list=who;
    opt=[];arg_list=[];
    pos=1;flag_thereismore=1;
    %  check if argument gives us a structure
    if flag_thereismore
        if isstruct(varargin{pos})
            opt=varargin{pos};
            pos=pos+1;
            flag_thereismore=nargin>(pos+1);
        end
    end
    %  check for name-variable pairs
    if flag_thereismore
        if ((nargin-pos-1)/2)==fix((nargin-pos-1)/2)
            arg_list=varargin(pos:end);
        else
            error('No of arguments is off.')
        end
    end  
    %  add option structure variables if they are part of the list
    if ~isempty(opt)
        for i1=1:numel(list)
            if isfield(opt, char(list{i1}))
                eval(horzcat(matlab.lang.makeValidName(char(list(i1))),'=opt.',char(list{i1}),';'));
            end
        end
    end
    %  add name-value pair arguments if they are part of the list
    if ~isempty(arg_list)
        for i1=1:numel(arg_list)/2
            if ismember(arg_list{(i1-1)*2+1},list) && isnumeric(arg_list{(i1-1)*2+2})
                eval(horzcat(arg_list{(i1-1)*2+1},'=',num2str(arg_list{(i1-1)*2+2}),';'));
            end
        end
    end 
end
%  PREP
n=numel(x);
m=1;                        %   no of output arguments
if isempty(D)
    D=NaN(MaxStepNo,MaxStepNo,n,m);
else
    siz=size(D);siz([1 2])=siz([1 2])+MaxStepNo;
    m=siz(4);
    D_new=NaN(siz);
    D_new(1:size(D,1),1:size(D,2),:,:)=D;
    D=D_new;
end
siz=size(D);
I=siz(1);
Q=inf(siz(1),n,m);
R=inf(siz(1),n,m);
J=NaN(m,n);
E=NaN(m,n);
S=ones(m,n);
%   size of increment before boundaries are considered
h=FinDiffRelStep.*max(abs(x),TypicalX);  
func_evals=0;
%  LOOP
for k=1:n   
    %  find final step size for given boundaries    
    if isfinite(lb(k)) || isfinite(ub(k))
        h(k)=min([h(k);abs(x(k)-lb(k))/2;abs(x(k)-ub(k))/2]);
        if h(k)==0
            error('evaluation on boundary not possible.')
        end
    end
    %  prep for this loop
    current_J=sum(~isnan(D(:,1,k,1)));
    for j=2:current_J
        Q(j,k,:)=abs(D(j,j,k,:)-D(j-1,j-1,k,:));
        R(j,k,:)=2*squeeze(Q(j,k,:))./squeeze(abs(D(j,j,k,:))+abs(D(j-1,j-1,k,:))+eps);
    end
    %  loop to increase depth
    flag_out=0;
    for j=(current_J+1):I
        if j==2
            if ~any(any(~(R(1:j-1,k,:)<toler),1),3)
                flag_out=1;
            end
        elseif j>2            
            if ~any(any(~(R(1:j-1,k,:)<toler),1) & all(Q(1:j-2,k,:)>=Q(2:j-1,k,:),1),3)
                flag_out=1;
            end
        end
        if flag_out==0
            %   only go on if one jacobian derivative element has not
            %   reduced its relative error below toler, and that element is
            %   still no suffering from numerical problems
            x_l=x;
            x_l(k)=x(k)-2^(-j)*h(k);
            x_r=x;
            x_r(k)=x(k)+2^(-j)*h(k);
            temp=(func(x_r)-func(x_l))/(2^(-j+1)*h(k));
            if m<numel(temp) && k==1 && j==1 % first evaluation, and we have more than one output argument
                m=numel(temp);
                Q(:,:,2:m)=inf;
                R(:,:,2:m)=inf;
            end
            D(j,1,k,1:m)=temp;%(func(x_r)-func(x_l))/(2^(-j+1)*h(k));
            func_evals=func_evals+2;
            for l=1:(j-1)
                D(j,l+1,k,:)=D(j,l,k,:)+(D(j,l,k,:)-D(j-1,l,k,:))/((4^l)-1);
            end
            if j>1
                Q(j,k,1:m)=abs(D(j,j,k,:)-D(j-1,j-1,k,:));
                R(j,k,1:m)=2*squeeze(Q(j,k,:))./squeeze(abs(D(j,j,k,:))+abs(D(j-1,j-1,k,:))+eps);
            end
        else
            break
        end
    end
    %  choose best fit for jacobian
    temp=Q(:,k,:);
    temp(~cumprod(Q(1:end-1,k,:)>=Q(2:end,k,:),1))=NaN;
    [err,pos]=min(temp,[],1);
    for i1=1:m
        J(i1,k)= D(pos(i1),pos(i1),k,i1);
        E(i1,k)= err(i1);
        S(i1,k)= pos(i1);
    end
end
end

%==================================================
% JACOBIAN LIMITED 
%==================================================
function  [J,h,func_evals,f_0,E,S,F,X] = jacobianlim(func,x,varargin)

%  SET OPTIONS
f_0=[];                     %   value at evaluation point (can be user-supplied
lb=-inf(size(x));           %   lower boundary of domain of X
ub=inf(size(x));            %   upper boundary of domain of X
FinDiffRelStep=eps^(1/3);   %   relative differentiation step size
TypicalX=1;                 %   sets lower limit of differentiation step 
%                               size, similar to the optimization routines 
%                               in matlab 
DecreaseStepSize=10;        %   factor that reduces the step size in limiting process
MaxStepNo=3;                %   max no of steps for the limiting process
toler=1e-6;                 %   absolute tolerances
%  DYNAMIC READ IN
if nargin>2    
    list=who;
    opt=[];arg_list=[];
    pos=1;flag_thereismore=1;
    %  check if argument gives us a structure
    if flag_thereismore
        if isstruct(varargin{pos})
            opt=varargin{pos};
            pos=pos+1;
            flag_thereismore=nargin>(pos+1);
        end
    end
    %  check for name-variable pairs
    if flag_thereismore
        if ((nargin-pos-1)/2)==fix((nargin-pos-1)/2)
            arg_list=varargin(pos:end);
        else
            error('No of arguments is off.')
        end
    end  
    %  add option structure variables if they are part of the list
    if ~isempty(opt)
        for i1=1:numel(list)
            if isfield(opt, char(list{i1}))
                eval(horzcat(matlab.lang.makeValidName(char(list(i1))),'=opt.',char(list{i1}),';'));
            end
        end
    end
    %  add name-value pair arguments if they are part of the list
    if ~isempty(arg_list)
        for i1=1:numel(arg_list)/2
            if ismember(arg_list{(i1-1)*2+1},list) && isnumeric(arg_list{(i1-1)*2+2})
                eval(horzcat(arg_list{(i1-1)*2+1},'=',num2str(arg_list{(i1-1)*2+2}),';'));
            end
        end
    end 
end
%  PREP
n=numel(x);
func_evals=0;
if isempty(f_0)
    f_0=func(x);
    func_evals=func_evals+1;
end
m=numel(f_0);
J=NaN(m,n);
E=NaN(m,n);
S=ones(m,n);
F=NaN(m,2*MaxStepNo+1,n);
X=NaN(n,2*MaxStepNo+1,n);
%   size of increment before boundaries are considered
h=FinDiffRelStep.*max(abs(x),TypicalX);  
%  ROUTINE FOR DIFFERENT PARAMETERS
for k=1:n
    %  find final step size for given boundaries    
    if isfinite(lb(k)) || isfinite(ub(k))
        h(k)=min([h(k);abs(x(k)-lb(k))/2;abs(x(k)-ub(k))/2]);
        if h(k)==0
            error('evaluation on boundary not possible.')
        end
    end
    %  prep for this loop
    D=NaN(m,MaxStepNo);
    D_r=NaN(m,MaxStepNo);
    D_l=NaN(m,MaxStepNo);
    Q=inf(m,MaxStepNo);
    R=zeros(m,MaxStepNo);
    Q_sided=inf(m,MaxStepNo);
    R_sided=zeros(m,MaxStepNo);
    %  first approximation of derivative
    [D(:,1),D_l(:,1),D_r(:,1),x_up,f_up,x_down,f_down]=deriv(func,f_0,x,h(k),k);
    func_evals=func_evals+2;
    F(:,[1 MaxStepNo+1 end],k)=[f_down(:) f_0(:) f_up(:)];
    X(:,[1 MaxStepNo+1 end],k)=[x_down(:) x(:) x_up(:)];
    Q_sided(:,1)=abs(D_l(:,1)-D_r(:,1));
    R_sided(:,1)=2*Q_sided(:,1).*(abs(D_l(:,1))+abs(D_r(:,1))+eps);
    %   see if first step is already sufficient
    i1=1;
    ind= find(R_sided(:,i1)<toler & ~isnan(D(:,i1)));
    J(ind,k)=D(ind,i1);
    E(ind,k)=Q_sided(ind,i1);
    S(ind,k)=i1;
    %  HIGHER LIMITS
    if any(isnan(J(:,k)))
        h(k)=h(k)/DecreaseStepSize;
        for i1=2:MaxStepNo
            %  evaluation
            [D(:,i1),D_l(:,i1),D_r(:,i1),x_up,f_up,x_down,f_down]=deriv(func,f_0,x,h(k),k);
            func_evals=func_evals+2;
            F(:,[i1 end-i1+1],k)=[f_down(:) f_up(:)];
            X(:,[i1 end-i1+1],k)=[x_down(:) x_up(:)];
            Q(:,i1)=abs(D(:,i1)-D(:,i1-1));
            R(:,i1)=2*Q(:,i1).*(abs(D(:,i1))+abs(D(:,i1-1))+eps);
            Q_sided(:,i1)=abs(D_l(:,i1)-D_r(:,i1));
            R_sided(:,i1)=2*Q_sided(:,i1).*(abs(D_l(:,i1))+abs(D_r(:,i1))+eps);
            %  current error tolerance is small enough
            ind= find( ~(Q(:,i1-1)<Q(:,i1)) & ( R(:,i1)<toler | R_sided(:,i1)<toler)  & ~isnan(D(:,i1)));%& isnan(J(:,order)) 
            J(ind,k)=D(ind,i1);
            E(ind,k)=Q(ind,i1);
            S(ind,k)=i1;
            %  error is growing and we take the last iterations value
            ind= find( ~(Q(:,i1-1)>=Q(:,i1)) & isnan(J(:,k)) & ~isnan(D(:,i1-1)));
            J(ind,k)=D(ind,i1-1);
            E(ind,k)=Q(ind,i1-1);
            S(ind,k)=i1-1;
            %  if error was growing, make sure future value growth reflects this
            Q(:,i1-1)=min(Q(:,i1-1),Q(:,i1));     
            if all(~isnan(J(:,k)))
                break
            else
                h(k)=h(k)/DecreaseStepSize;
            end
        end
    end
    %   collecting results
    ind=find(isnan(J(:,k)));
    J(ind,k)=D(ind,end);
    E(ind,k)=Q(ind,end);
    S(ind,k)=MaxStepNo;    
end
end
function [out,out_l,out_r,x_up,f_up,x_down,f_down]=deriv(func,f_0,x,h,order)
%   right point
x_up=x;
x_up(order)=x(order)+h;
%   right evaluation
f_up=func(x_up);f_up=f_up(:);
%   left point
x_down=x;
x_down(order)=x(order)-h;
%   left evaluation
f_down=func(x_down);f_down=f_down(:);
%   central derivative
out=diff([f_down(:) f_up]')'/(2*h);
%   left derivative
out_l=diff([f_down f_0(:)]')'/h;
%   right derivative
out_r=diff([f_0(:) f_up]')'/h;
end

