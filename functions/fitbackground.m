%
% FITBACKGROUND Fit the background function in a signal
%
%   [B,lambda,param,tstart] = FITBACKGROUND(V,t,@model)
%   Fits the background (B) and the modulation depth (lambda) to a
%   time-domain signal (V) and time-axis (t) based on a given time-domain
%   parametric model (@model). The fitted parameters of the model are
%   returned as a last output argument. The optimal fitting start time (tstart)
%   is computed automatically using the backgroundstart() function.
%
%   [B,lambda,param] = FITBACKGROUND(V,t,@model,tstart)
%   The time at which the background starts to be fitted can be passed as a
%   an additional argument.
%
%   [B,lambda,param] = FITBACKGROUND(V,t,@model,[tstart tend])
%   The start and end times of the fitting can be specified by passing a
%   two-element array as the argument. If tend is not specified, the end of
%   the signal is selected as the default.
%
%   [B,lambda,param] = FITBACKGROUND(...,'Property',Value)
%   Additional (optional) arguments can be passed as property-value pairs.
%
% The properties to be passed as options can be set in any order.
%
%   'ModDepth' - Fixes the modulation depth to a user-defined value instead
%                of fitting it along the background.
%   'LogFit' - Specifies whether to fit the log of the signal (default: false)
%   'InitialGuess' - Array of initial values for the fit parameters
%   'Solver' - Optimization solver used for the fitting
%             'lsqnonlin' - Non-linear constrained least-squares (toolbox)
%             'nlsqbnd'   - Non-linear constrained least-squares (free)
%

% This file is a part of DeerLab. License is MIT (see LICENSE.md).
% Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

function [B,ModDepth,FitParam,tFitRange] = fitbackground(V,t,bgmodel,tFitRange,varargin)


if nargin<3
    error('Not enough input arguments. At least three are needed: V, t, and background model.');
end

if nargin<4
    tstart = backgroundstart(V,t,bgmodel);
    tend = t(end);
    tFitRange = [tstart tend];
elseif ischar(tFitRange)
    varargin = [{tFitRange} varargin];
    tstart = backgroundstart(V,t,bgmodel);
    tend = t(end);
    tFitRange = [tstart tend];
    
elseif length(tFitRange) == 1
    tFitRange(2) = t(end);
elseif length(tFitRange) > 2
    error('The 4th argument cannot exceed two elements.')
end

if tFitRange(2)<tFitRange(1)
    error('The fit start time cannot exceed the fit end time.')
end

if ~isa(bgmodel,'function_handle')
    error('The background model (third input) must be a valid function handle.')
end

% Parse optional inputs
[LogFit,InitialGuess,ModDepth,Solver] = parseoptional({'LogFit','InitialGuess','ModDepth','Solver'},varargin);

if isempty(LogFit)
    LogFit = false;
end
if ~isempty(ModDepth)
    if ~isnumeric(ModDepth) || length(ModDepth)>1
        error('Modulation depth in ''ModDepth'' must be a scalar.')
    end
    if ModDepth>1 || ModDepth<0
        error('Fixed modulation depth in ''ModDepth'' must be in the range 0 to 1.')
    end
end
fitModDepth = isempty(ModDepth);

OptimizationToolboxInstalled = optimtoolbox_installed();

if isempty(Solver)
    if OptimizationToolboxInstalled 
        Solver = 'lsqnonlin';
    else
        Solver = 'fminsearchcon';
    end
else
    validateattributes(Solver,{'char'},{'nonempty'},mfilename,'Solver')
end

SolverList = {'lsqnonlin','nlsqbnd','fminsearchcon'};
validatestring(Solver,SolverList);

if strcmp(Solver,'nlsqbnd') && ~ispc
   error('The ''nlsqbnd'' solver is only available for Windows systems.') 
end

if strcmp(Solver,'lsqnonlin') && ~OptimizationToolboxInstalled 
    error('The ''lsqnonlin'' solver requires the Optimization Toolbox.')
end

% Ensure real part is used
V = real(V);
% Validate inputs
validateattributes(InitialGuess,{'numeric'},{'2d'},mfilename,'InitialGuess')
validateattributes(tFitRange,{'numeric'},{'2d','nonempty'},mfilename,'FitDelimiter')
validateattributes(V,{'numeric'},{'2d','nonempty'},mfilename,'Data')
validateattributes(t,{'numeric'},{'2d','nonempty','increasing'},mfilename,'t')

% Use column vectors
t = t(:);
V = V(:);

%--------------------------------------------------------------------------

% Find the position to limit fit
FitStartTime = tFitRange(1);
FitEndTime = tFitRange(2);
[~,FitStartPos] = min(abs(t - FitStartTime));
[~,FitEndPos] = min(abs(t - FitEndTime));

% Limit the time axis and the data to fit
tfit = t(FitStartPos:FitEndPos);
FitData = V(FitStartPos:FitEndPos);

% Construct objective function for minimization
% Fit signal or log(signal)
Fgmodel = @(p,lambda)(1 - lambda + eps)*bgmodel(tfit,p);
if LogFit
    residuals = @(p,lambda) sqrt(1/2)*(log(Fgmodel(p,lambda)) - log(FitData));
else
    residuals = @(p,lambda) sqrt(1/2)*(Fgmodel(p,lambda) - FitData);
end
if fitModDepth
    ObjFcnVec = @(param) residuals(param(1:end-1),param(end));
    ObjFcn = @(param) norm(residuals(param(1:end-1),param(end)))^2;
else
    ObjFcnVec = @(param) residuals(param,ModDepth);
    ObjFcn = @(param) norm(residuals(param,ModDepth))^2;
end

% Initialize bounds and initial parameter values
info = bgmodel();
parinfo = info.parameters;
Ranges = cat(1,parinfo.range);
lowerBounds = Ranges(:,1);
upperBounds = Ranges(:,2);
if fitModDepth
    lowerBounds(end+1) = 0;
    upperBounds(end+1) = 1;
end
if ~isempty(InitialGuess)
    StartParameters = InitialGuess;
else
    StartParameters =  [parinfo.default];
    if fitModDepth
        StartParameters(end+1) = 0.5;
    end
end

% Solve the constrained nonlinear minimization problem
switch lower(Solver)
    case 'lsqnonlin'
        % Prepare minimization problem solver
        solveropts = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off',...
            'MaxIter',8000,'MaxFunEvals',8000,...
            'TolFun',1e-10,'DiffMinChange',1e-8,'DiffMaxChange',0.1);
        % Run solver
        FitParam = lsqnonlin(ObjFcnVec,StartParameters,lowerBounds,upperBounds,solveropts);
    case 'nlsqbnd'
        % Prepare minimization problem solver
        solveropts = optimset('Algorithm','trust-region-reflective','Display','off',...
            'MaxIter',8000,'MaxFunEvals',8000,...
            'TolFun',1e-20,'TolCon',1e-20,...
            'DiffMinChange',1e-8,'DiffMaxChange',0.1);
        % Run solver
        FitParam = nlsqbnd(ObjFcnVec,StartParameters,lowerBounds,upperBounds,solveropts);
        % nlsqbnd returns a column, transpose to adapt to row-style of MATLAB solvers
        FitParam = FitParam.';
    case 'fminsearchcon'
        % Prepare minimization problem solver
        solverOpts = optimset('Algorithm','trust-region-reflective','Display','off',...
            'MaxIter',8000,'MaxFunEvals',8000,...
            'TolFun',1e-20,'TolCon',1e-20,...
            'DiffMinChange',1e-8,'DiffMaxChange',0.1);
        % Run solver
        FitParam = fminsearchcon(ObjFcn,StartParameters,lowerBounds,upperBounds,[],[],[],solverOpts);
end

% Extract the fitted modulation depth
if fitModDepth
    ModDepth = FitParam(end);
    FitParam = FitParam(1:end-1);
end

% Extrapolate fitted background to whole time axis
B = bgmodel(t,FitParam);

% Ensure data is real
B = real(B);
B = B(:);

end
