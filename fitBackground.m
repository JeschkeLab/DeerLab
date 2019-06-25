% fitBackground - Background fitting function of DeerAnalysis
% handles.background
% 0  fractal, n variable,  exp(-k*t^(n/3))
% 1  n-dimensional, n fixed, exp(-k*t^(n/3))
% 2  three-dimensional, exp(-k*t)
% 3  polynomial
% 4  user-defined function in handles.bckg_fct or numerical background in handles.bckg_data
% 5  homogeneous background
%
% texp      full time range
% t_fit     time range for background fit
% data_fit  time-domain data over t_fit

function [Background,FitResults]=fitBackground(FitData,TimeAxis,FitTimeAxis,BckgModel,ModelParam)

if nargin<3
  error('Not enough input arguments.')
end

if nargin<4
  BckgModel = 'exponential';
end

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
    
  case 4
    bckg0=exp(polyval(handles.polynomial,abs(FitTimeAxis)));
    LinearLogFit=[1 0.8];
    fminResults=fminsearch(@rms_ubckg,LinearLogFit,[],FitData,bckg0);
    bckg1=exp(polyval(handles.polynomial,abs(TimeAxis)));
    Background=fminResults(2)*exp(fminResults(1)*log(bckg1));
    density=fminResults(1)/handles.calib_density;
    handles.bckg_dens=density;
    
  case 5
    fminResults(1)=handles.man_k;
    fminResults(2)=1;
    Background=decaynD(fminResults,TimeAxis,handles.hom_dim);
    Background=(1-handles.man_depth)*Background;
    % bckg=(1-handles.man_depth)*exp(-dens*targ);
    density=handles.man_k;
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

