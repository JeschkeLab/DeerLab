function [Background,FitResults] = fitBackground(FitData,TimeAxis,FitTimeAxis,BckgModel,ModelParam)

if nargin<3
    error('Not enough input arguments.')
end

if nargin<4 || isempty(BckgModel)
    BckgModel = 'exponential';
elseif isa(BckgModel,'fittype')
    CustomFitModel = BckgModel;
    BckgModel = 'custom';
elseif ~isa(BckgModel,'fittype') && ~isa(BckgModel,'char')
    error('BckgModel must be a valid ''fittype'' or ''char'' class variable.')
else
    validateattributes(BckgModel,{'char'},{'nonempty'},mfilename,'BckgModel')
    allowedInput = {'fractal','exponential','polyexp','polynomial'};
    validatestring(BckgModel,allowedInput);
end
if iscolumn(TimeAxis)
    TimeAxis = TimeAxis';    
end
if iscolumn(FitData)
    FitData = FitData';
    DataIsColumn = true;
else
    DataIsColumn = false;
end
if iscolumn(FitTimeAxis)
    FitTimeAxis = FitTimeAxis';
end
validateattributes(FitData,{'numeric'},{'2d','nonempty'},mfilename,'FitData')
validateattributes(TimeAxis,{'numeric'},{'2d','nonempty','nonnegative','increasing'},mfilename,'TimeAxis')


%Start with a linear fit of log-data
PolynomialFit=polyfit(FitTimeAxis,log(FitData),1);
LinearLogFit=[-PolynomialFit(1) 1];

switch BckgModel
    
    case 'fractal'
        %Use linear log-fit as initial gues
        LinearLogFit = [-PolynomialFit(1) 1 3];
        %Fit amplitude, decay rate and fractal dimension of stretched exp.
        fminResults = fminsearch(@minimizeStretchExp,LinearLogFit,[],FitTimeAxis,FitData);
        %Limit resolution to fourth digit to avoid precission errors
        fminResults = round(fminResults,4);
        %Compute stretched exponential background
        Background = fminResults(2)*exp(-abs(fminResults(1)*TimeAxis).^(fminResults(3)/3));
        FitResults.DecayRate = fminResults(1);
        FitResults.FractalDimension = fminResults(3);
        
    case 'exponential'
        %Compute exponential background from linear log-fit
        Background = exp(polyval(PolynomialFit,abs(TimeAxis)));
        FitResults.DecayRate = LinearLogFit(1);
        
    case 'polynomial'
        %Fit polynomial function
        PolynomialFit = polyfit(FitTimeAxis,FitData,ModelParam);
        %Compute polynomial background
        Background = polyval(PolynomialFit,abs(TimeAxis));
        FitResults.polynomial = PolynomialFit;
        
    case 'polyexp'
        %Fit polynomial function to log-data
        PolynomialFit = polyfit(FitTimeAxis,log(FitData),1);
        %Compute polynomial background
        Background = exp(polyval(PolynomialFit,abs(TimeAxis)));
        FitResults.polynomial = PolynomialFit;
    case 'custom'
        %Fit using the user-provided model
        FitObject = fit(FitTimeAxis',FitData',CustomFitModel);
        %Compute fitted background
        Background = feval(FitObject,TimeAxis)';
end

%Ensure data is real
Background=real(Background);
if DataIsColumn
   Background = Background'; 
end
return

%Cost function for fractal stretched exponential fit
function RMSD = minimizeStretchExp(StretchedExpParam,TimeAxis,Obervation)
if StretchedExpParam(1)<0
    RMSD=1.0e10;
    return;
end
Prediction = StretchedExpParam(2)*exp(-abs(StretchedExpParam(1)*TimeAxis).^(StretchedExpParam(3)/3));
RMSD = 1/2*norm(Prediction - Obervation)^2;
return

