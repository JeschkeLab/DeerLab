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
%   'LogFit' - Specifies whether to fit the log of the signal (default: false)
%
%   'InitialGuess' - Array of initial values for the fit parameters
%

% This file is a part of DeerAnalysis. License is MIT (see LICENSE.md). 
% Copyright(c) 2019: Luis Fabregas, Stefan Stoll, Gunnar Jeschke and other contributors.

function [B,ModDepth,FitParam,FitDelimiter] = fitbackground(Data,t,BckgModel,FitDelimiter,varargin)


if ~license('test','optimization_toolbox')
   error('DeerAnaysis could not find a valid licence for the Optimization Toolbox. Please install the add-on to use fitbackground.')
end

if nargin<3
    error('Not enough input arguments. At least three are needed: V, t, and background model.');
end

if nargin<4
    tstart = backgroundstart(Data,t,BckgModel);
    tend = t(end);
    FitDelimiter = [tstart tend];
elseif ischar(FitDelimiter)
    varargin = [{FitDelimiter} varargin];
    tstart = backgroundstart(Data,t,BckgModel);
    tend = t(end);
    FitDelimiter = [tstart tend]; 
    
elseif length(FitDelimiter) == 1
    FitDelimiter(2) = t(end);
elseif length(FitDelimiter) > 2
    error('The 4th argument cannot exceed two elements.')
end

if FitDelimiter(2)<FitDelimiter(1)
    error('The fit start time cannot exceed the fit end time.')
end

tstart = FitDelimiter(1);

if ~isa(BckgModel,'function_handle')
    error('The background model must be a valid function handle.')
end

if ~iscolumn(t)
    t = t.';
end

%Parse optiona inputs
[LogFit,InitialGuess,ModDepth] = parseoptional({'LogFit','InitialGuess','ModDepth'},varargin);

if isempty(LogFit)
    LogFit = false;
end
if length(ModDepth)>1
    error('FixLambda must be a scalar.')
end
if ~isempty(ModDepth)
    if ModDepth>1 || ModDepth<0
        error('Fixed modulation depth must be in the range 0 to 1.')
    end
end
DataIsColumn = iscolumn(Data);
if ~DataIsColumn
    Data = Data.';
end
Data = real(Data);
validateattributes(InitialGuess,{'numeric'},{'2d'},mfilename,'InitialGuess')
validateattributes(FitDelimiter,{'numeric'},{'2d','nonempty'},mfilename,'FitDelimiter')
validateattributes(Data,{'numeric'},{'2d','nonempty'},mfilename,'Data')
validateattributes(t,{'numeric'},{'2d','nonempty','increasing'},mfilename,'t')


%Convert time step to microseconds if given in nanoseconds
if isnanosecond(t)
    t = t/1000; % ns->us
    FitDelimiter = FitDelimiter/1000; % ns->us
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%Find the position to limit fit
FitStartTime = FitDelimiter(1);
FitEndTime = FitDelimiter(2);
[~,FitStartPos] = min(abs(t - FitStartTime));
[~,FitEndPos] = min(abs(t - FitEndTime));

%Limit the time axis and the data to fit
Fitt = t(FitStartPos:FitEndPos);
FitData = Data(FitStartPos:FitEndPos);

%Use absolute time scale to ensure proper fitting of negative-time data
Fitt = abs(Fitt);

%Prepare minimization problem solver
solveropts = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective','Display','off',...
    'MaxIter',8000,'MaxFunEvals',8000,...
    'TolFun',1e-10,'DiffMinChange',1e-8,'DiffMaxChange',0.1);

%Construct cost functional for minimization
if LogFit
    %Fit in log-space
    if isempty(ModDepth)
        %Fit modulation depth...
        CostFcn = @(param)(sqrt(1/2)*(log((1 - param(1) + eps)*BckgModel(Fitt,param(2:end))) - log(FitData)));
    else
        %Use user-given modulation depth...
        CostFcn = @(param)(sqrt(1/2)*(log((1 - ModDepth + eps)*BckgModel(Fitt,param(1:end))) - log(FitData)));
    end
else
    %Fit in linear space
    if isempty(ModDepth)
        %Fit modulation depth...
        CostFcn = @(param)(sqrt(1/2)*((1 - param(1))*BckgModel(Fitt,param(2:end)) - FitData));
    else
        %Use user-given modulation depth...
        CostFcn = @(param)(sqrt(1/2)*((1 - ModDepth)*BckgModel(Fitt,param(1:end)) - FitData));
    end
end

if isempty(ModDepth)
    %Initiallize StartParameters (1st element is modulation depth)
    LowerBounds(1) = 0;
    UpperBounds(1) = 1;
end
pos = 1 - length(ModDepth);
%Get information about the time-domain parametric model
Info = BckgModel();
Ranges =  [Info.parameters(:).range];
LowerBounds(1+pos:pos + Info.nparam) = Ranges(1:2:end-1);
UpperBounds(1+pos:pos + Info.nparam) = Ranges(2:2:end);
if ~isempty(InitialGuess)
    StartParameters = InitialGuess;
else
    if isempty(ModDepth)
        StartParameters(1) = 0.5;
    end
    StartParameters(1+pos:pos + Info.nparam) =  [Info.parameters(:).default];
end

FitParam = lsqnonlin(CostFcn,StartParameters,LowerBounds,UpperBounds,solveropts);

if isempty(ModDepth)
    %Get the fitted modulation depth
    ModDepth = FitParam(1);
end
%Remove the modulation depth from the fit parameters
FitParam = FitParam(1+pos:end);

%Extrapolate fitted background to whole time axis
B = BckgModel(abs(t),FitParam);

%Ensure data is real
B = real(B);
if ~DataIsColumn
    B = B';
end

end
