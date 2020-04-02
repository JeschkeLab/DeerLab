function [Pfit,Vfit,Bfit] = fitsignal(V,t,r,dd_model,bg_model,ex_model)

% Use simple default models if not specified
if nargin<5
    ex_model = @ex_4pdeer;
end
if nargin<4
    bg_model = @bg_strexp;
end

if isa(dd_model,'function_handle')
    % Get information about distance distribution model
    info = dd_model();
    modelidx = ones(1,info.nparam)*1;
    param0 = [info.parameters.default];
    range = [info.parameters.range];
    lower = range(1:2:end-1);
    upper = range(2:2:end);
else
    param0 = [];
    lower = [];
    upper = [];
    modelidx = [];
end

% Get information about background model
info = bg_model();
modelidx = [modelidx ones(1,info.nparam)*2];
param0 = [param0 info.parameters.default];
range = [info.parameters.range];
lower = [lower range(1:2:end-1)];
upper = [upper range(2:2:end)];

% Get information about experiment model
info = ex_model(t);
modelidx = [modelidx ones(1,info.nparam)*3];
param0 = [param0 info.parameters.default];
range = [info.parameters.range];
lower = [lower range(1:2:end-1)];
upper = [upper range(2:2:end)];

% Call fitparamodel to fit the parameterset
parfit = fitparamodel(V,@Vmodel,t,param0,'Lower',lower,'Upper',upper);

% Get the model fits
[Vfit,Bfit,Pfit] = Vmodel(t,parfit);


    % Function for our general multipulse DEER signal model
    function [Vfit,Bfit,Pfit] = Vmodel(t,par)
        
        % Extreact parameter subsets
        theta_dd = par(modelidx == 1);
        theta_bg = par(modelidx == 2);
        theta_ex = par(modelidx == 3);

        % Prepare the background model
        Bbasis = @(t)bg_model(t,theta_bg);

        % Get the dipolar kernel and background according to experiment
        [K,Bfit] = ex_model(t,r,theta_ex,Bbasis);
        
        % Get the distance distribution from the model or via regularization
        if isa(dd_model,'function_handle')
            % Parametric distribution
            Pfit = dd_model(r,theta_dd);
        else
            % Parameter-free distribution
            Pfit = fitregmodel(V,K,r,'tikh','aic');
        end
        
        % Finally, compute the dipolar signal
        Vfit = K*Pfit;
        
    end

    

end
