function [pass,maxerr] = test(opt)

% Test selregparam with Tikhonov regularization using golden search

t = linspace(0,3,200);
r = linspace(2,6,100);
P = rd_onegaussian(r,[3,0.5]);
K = dipolarkernel(t,r);
S = K*P;

[alphaopt,fcns,alphas] = selregparam(S,K,r,'tikhonov',{'aic','gml','gcv'},'Search','grid');

% Pass: similar regularization parameter values are found
pass = all(abs(diff(alphaopt)) < 1e-2);
 
maxerr = max(diff(alphaopt));

if opt.Display

    subplot(131)
    plot(log(alphas),log(fcns{1}+ 1e5),'.')
    xlabel('log_{10}(\alpha)')
    yabel('AIC')
    grid on, axis tight, box on
    
    subplot(132)
    plot(log(alphas),log(fcns{2}),'.')
    xlabel('log_{10}(\alpha)')
    yabel('GML')
    grid on, axis tight, box on
    
    subplot(133)
    plot(log(alphas),log(fcns{3}),'.')
    xlabel('log_{10}(\alpha)')
    yabel('GCV')
    grid on, axis tight, box on
    
end


end