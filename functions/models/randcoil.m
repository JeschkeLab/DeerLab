function output = randcoil(r,param)

% Random-coil model for an unfolded peptide/protein, end-to-end distance
% distribution is approximated by a Gaussian coil with proper mean
% distance, which is good for sufficiently large N
% N. C. Fitzkee, G. D. Rose, PNAS 2004, 101(34), 12497-12502
%
% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  N      50      2              1000    number of residues between labels, including labeled residues
% par(2)  nu     0.602   0.33           1       scaling exponent

R0=0.198; % 1.98 Å per residue
nParam = 2;

if nargin==0
    %If no inputs given, return info about the parametric model
    info.Model  = 'Random-coil model';
    info.Equation  = ['3/(2*pi*3/(2*pi*6.r0*N*',char(957),')^(3/2)))^(3/2)*4*pi*r.^2*exp(-3*r.^2/(3/(2*pi*6.r0*N*',char(957),')^(3/2))'];
    info.nParam  = nParam;
    info.parameters(1).name = 'Chain members N';
    info.parameters(1).range = [2 1000];
    info.parameters(1).default = 50;
    info.parameters(1).units = '';
    
    info.parameters(2).name = ['Scaling exponent',char(957)];
    info.parameters(2).range = [0.33 1];
    info.parameters(2).default = 0.602;
    info.parameters(2).units = '';
    
    output = info;
    
elseif nargin == 2
    
    %If user passes them, check that the number of parameters matches the model
    if length(param)~=nParam
        error('The number of input parameters does not match the number of model parameters.')
    end
    
    SquareDist = 6*R0*param(1)^param(2)^2; %mean square end-to-end distance from radius of gyration
    normFact = 3/(2*pi*SquareDist)^(3/2); % normalization prefactor
    ShellSurf = 4*pi*r.^2; % spherical shell surface
    Gaussian = exp(-3*r.^2/(2*SquareDist));
    Distribution = normFact*ShellSurf.*Gaussian;
    
    %Normalize integral
    Distribution = Distribution/sum(Distribution)/mean(diff(r));
    output = Distribution;
    
else
    
    %Else, the user has given wrong number of inputs
    error('Model requires two input arguments.')
end

return

