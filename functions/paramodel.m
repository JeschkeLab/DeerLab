%
% PARAMODEL Converts a function handle to a DeerAnalysis parametric model
%
%       g = PARAMODEL(f,n)
%       Converts the input function handle (f) to a valid DeerAnalysis
%       parametric model compatible with the fit functions. The number of
%       parameters in the resulting parametric model must be specified by
%       (n). The initial guess values of all parameters are set to zero.
%       The resulting parametric model g is returned in as a function
%       handle.
%
%       g = PARAMODEL(f,param0)
%       If an array of inital guess values are passed, these are set into
%       the resulting parametric model. The number of parameters is computed
%       from the length of the array.
%
% Copyright(C) 2019  Luis Fabregas, DeerAnalysis2
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License 3.0 as published by
% the Free Software Foundation.


function model = paramodel(handle,param0,lower,upper)


%Check which king of input is being use
if length(param0)==1
    nParam = param0;
    param0 = zeros(nParam,1);
else
    nParam = length(param0);
end

if nargin<3
    lower = zeros(nParam,1) - realmax;
    upper = zeros(nParam,1) + realmax;
end
if nargin<4
    upper = zeros(nParam,1) + realmax;
end

if any(length(param0)~=length(upper) & length(param0)~=length(lower))
   error('The Inital guess and upper/lower boundaries must have equal length); ') 
end

model = @myparametricmodel;

%Define the raw structure of the DeerAnalysis parametric model functions
    function Output = myparametricmodel(t,param)
        
        if nargin==0
            %If no inputs given, return info about the parametric model
            info.Model  = 'Automatically converted parametric model';
            info.Equation  = 'custom';
            info.nParam  = nParam;
            for i=1:nParam
                info.parameters(i).name = ' ';
                info.parameters(i).range = [lower(i) upper(i)];
                info.parameters(i).default = param0(i);
                info.parameters(i).units = ' ';
            end
            
            Output = info;
            
        elseif nargin == 2
            
            %If user passes them, check that the number of parameters matches the model
            if length(param)~=nParam
                error('The number of input parameters does not match the number of model parameters.')
            end
            
            %If necessary inputs given, compute the model distance distribution
            t = abs(t);
            Output = handle(t,param);
            if ~iscolumn(Output)
                Output = Output';
            end
        else
            %Else, the user has given wrong number of inputs
            error('Model requires two input arguments.')
        end
        
    end

end
