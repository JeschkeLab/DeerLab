% fitBackground - Background fitting function of DeerAnalysis
%
% Usage:
%   Fits the FitData with corresponding FitTimeAxis to a certain background
%   model. The background is  then extrapolated to the input TimeAxis.
%
%   fitBackground(FitData,TimeAxis,FitTimeAxis,BckgModel,ModelParam)
%   
%   BckgModel   'fractal'       stretched exponential 
%                               B(t) = exp(-k*t^(d/3))
%   
%               'exponential'   exponential function
%                               B(t) = exp(-k*t)
%   
%               'polynomial'    N-th order polynomial function 
%                               B(t) = sum_i^N  c_i*t^i 
%                               To give as parameter:
%                               ModelParam = N;
%   
%               'polyexp'       exponential function fitted as log
%                               log(B(t)) = -k*t
%
%    If not specified, BckgModel will be set to 'exponential'.
% 
% Luis Fabregas 2019, DeerAnalysis 
%
function [Background,FitResults]=fitBackground(FitData,TimeAxis,FitTimeAxis,BckgModel,ModelParam,Averaging)

if nargin<3
  error('Not enough input arguments.')
end

if nargin<4
  BckgModel = 'exponential';
end

if nargin<6
  Averaging = false;
end

%Start with a linear fit of log-data
PolynomialFit=polyfit(FitTimeAxis,log(FitData),1); 
LinearLogFit=[-PolynomialFit(1) 1];

if Averaging
Prediction = polyval(PolynomialFit,FitTimeAxis);
clear RMSD
for j=1:75
AverFitFitData = log(FitData);
NAverages = j;
for i=1:NAverages
AverFitFitData = movmean(AverFitFitData,7);
end
RMSD(j) = 1/2*norm(Prediction - AverFitFitData)^2;
end
[~,NAveragesOpt]  = min(RMSD);

AverFitFitData = FitData;
% NAveragesOpt
for i=1:NAveragesOpt
AverFitFitData = movmean(AverFitFitData,7);
end
FitData = AverFitFitData;

%Start with a linear fit of log-data
PolynomialFit=polyfit(FitTimeAxis,log(FitData),1); 
LinearLogFit=[-PolynomialFit(1) 1];
end


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

end

%Ensure data is real
Background=real(Background);

return

%Cost function for fractal stretched exponential fit  
function RMSD = minimizeStretchExp(StretchedExpParam,TimeAxis,Obervation)
if StretchedExpParam(1)<0
  RMSD=1.0e10;
  return;
end;
Prediction = StretchedExpParam(2)*exp(-abs(StretchedExpParam(1)*TimeAxis).^(StretchedExpParam(3)/3));
RMSD = 1/2*norm(Prediction - Obervation)^2;
return

