function bg_exvol_alpha()

% Precalculates DEER background reduction factor alpha(d)
% See
%    Kattnig et al
%    J.Phys. Chem. B, 117, 16542 (2013)
%    https://doi.org/10.1021/jp408338q
% The background reduction factor alpha(d) is defined in Eq.(18)
% For large d, one can use the limiting expression
%    alpha = (3/2/pi)*(2*pi/3-sqrt(3)./d); 
% as an excellent approximation (error at d

filename = 'bg_exvol.mat';
saveData = false;
doPlotting = true;

% Set up dR range
%-------------------------------------------------------------------------------
% dR = A*t/R^3, where t is time, R is excluded-volume radius, and A is the
% dipolar constant (in units compatible with t and R)
dRlin = 0:0.05:20;
dRlog = 10.^(1:0.05:3);
dRlog(dRlog<max(dRlin)) = [];
dR = [dRlin dRlog];

% Evaluate reduction factor alpha as a function of dR
%-------------------------------------------------------------------------------
h_ = h(dR);
K_ = KK(dR);
alpha = (3/2/pi)*(h_-sqrt(3)./dR.*K_);
alpha(dR==0) = 0;

% Plotting
%-------------------------------------------------------------------------------
% Reproduces Fig.1 from Kattnig et al.
if doPlotting
    alphalim = (3/2/pi)*(2*pi/3-sqrt(3)./dR);
    
    subplot(2,2,1)
    plot(dR,alpha,'.-',dR,alphalim);
    xlabel('dR');
    ylabel('\alpha');
    grid on
    %xlim([0 40]);
  
    % Plot entire dR range
    subplot(2,2,2);
    plot(dR,alpha-alphalim);
    xlabel('dR');
    ylabel('1-\alpha');
    grid on
    
    subplot(2,2,3);
    plot(dR,h_,'.-',dR,2*pi/3*ones(size(dR)));
    xlabel('dR');
    ylabel('h');
    
    subplot(2,2,4);
    plot(dR,K_,'.-',dR,ones(size(dR)));
    xlabel('dR');
    ylabel('K');
end


% Save data
%-------------------------------------------------------------------------------
if saveData
    exvol.dR = dR;
    exvol.alpha = alpha;
    save(filename,'exvol');
end

end

function y = KK(d)
q = sqrt(6*d/pi);
y = 1 - (cos(d).*fresnelC(q)+sin(d).*fresnelS(q))./q;
y(y==0) = 0;
end

function y = h(d)
y = zeros(size(d));
for k = 1:numel(d)
  fprintf('%d/%d - %g\n',k,numel(d),d(k));
  y(k) = integral(@(x)(1-x.^2).*Si((1-x.^2)*d(k)),0,sqrt(3));
end
end

function y = Si(t)
y = zeros(size(t));
for k = 1:numel(t)
  y(k) = integral(@(x)sin(x)./x,0,t(k));
end
y(y==0) = 0;
end
