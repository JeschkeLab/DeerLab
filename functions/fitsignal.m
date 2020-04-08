function [Pfit,Vfit,Bfit,parfit] = fitsignal(Vexp,t,r,dd_model,bg_model,ex_model)

if nargin<5
    error('At least five inputs (V,t,r,@dd,@bg) need to be specified.');
end

% Use simple default experiment model if not specified
if nargin<6
    ex_model = @exp_4pdeer;
end

% Get information about distance distribution parameters
if isa(dd_model,'function_handle')
    [par0_dd,lower_dd,upper_dd] = getmodelparams(dd_model);
elseif isempty(dd_model)
    par0_dd = [];
    lower_dd = [];
    upper_dd = [];
else
    error('Distribution model (4th input) must either be a function handle (for a parametric model) or [] (for a parameter-free distribution).')
end
Ndd = numel(par0_dd);

% Get information about background parameters
if isa(bg_model,'function_handle')
    [par0_bg,lower_bg,upper_bg] = getmodelparams(bg_model);
elseif isempty(bg_model)
    par0_bg = [];
    lower_bg = [];
    upper_bg = [];
else
    error('Background model (5th input) must either be a function handle, or [] if no background should be fitted.')
end
Nbg = numel(par0_bg);

% Get information about experiment parameters
if isa(ex_model,'function_handle')
    [par0_ex,lower_ex,upper_ex] = getmodelparams(ex_model);
else
    error('Experiment model (6th input) must be a function handle.');
end
Nex = numel(par0_ex);

% Combine all parameters into a single vector
par0 = [par0_dd par0_bg par0_ex];
lower = [lower_dd lower_bg lower_ex];
upper = [upper_dd upper_bg upper_ex];
modelidx = [ones(1,Ndd) ones(1,Nbg)*2 ones(1,Nex)*3];
ddidx = modelidx==1;
bgidx = modelidx==2;
exidx = modelidx==3;

% Fit the parameters
parfit_ = fitparamodel(Vexp,@Vmodel,t,par0,'Lower',lower,'Upper',upper);

% Calculate the fitted signal, background, and distribution
[Vfit,Bfit,Pfit] = Vmodel(t,parfit_);

% Return fitted parameter in structure
parfit.dd = parfit_(ddidx);
parfit.bg = parfit_(bgidx);
parfit.ex = parfit_(exidx);


    % General multi-pathway DEER signal model function
    function [V,B,P] = Vmodel(t,par)
        
        % Extract parameter subsets
        par_dd = par(ddidx);
        par_bg = par(bgidx);
        par_ex = par(exidx);
        
        % Calculate the experiment kernel and the background
        if isa(bg_model,'function_handle')
            Bfcn = @(t) bg_model(t,par_bg);
        else
            Bfcn = [];
        end
        [K,B] = ex_model(t,r,par_ex,Bfcn);
        
        % Get the distance distribution from the model or via regularization
        if isa(dd_model,'function_handle')
            P = dd_model(r,par_dd);
        else
            P = fitregmodel(Vexp,K,r,'tikh','aic');
        end
        
        % Compute the dipolar signal
        V = K*P;
        
    end
    
end

function [par0,lo,up] = getmodelparams(model)

info = model();
par0 = [info.parameters.default];
range = [info.parameters.range];
lo = range(1:2:end-1);
up = range(2:2:end);

end
