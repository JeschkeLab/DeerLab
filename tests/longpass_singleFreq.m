function [err,data,maxerr] = test(opt,olddata)

%=======================================
% Check TV regularization
%=======================================
Dimension = 200;
dt = 0.008;
t = linspace(0,dt*Dimension,Dimension);
r = time2dist(t);
P = rd_onegaussian(r,[3,0.5]);

K = dipolarkernel(t,r);
DipEvoFcn = K*P;


rESEEM = 1.5;
EseemFreq = 52.04/(rESEEM^3);
ESEEM = 0.5*(2 + exp(-7*t).*cos(2*pi*EseemFreq*t));
S = DipEvoFcn.*ESEEM';

Filtered = longpass(t,S,1.5);

L = regoperator(Dimension,2);
RegParam = regparamrange(K,L);
RegParam2 = selregparam(Filtered,K,L,'tikhonov','gml');
Result = fitregmodel(Filtered,K,r,L,'tikhonov',RegParam2,'Solver','fnnls');

error = abs(Result - P);
err(1) = any(error>3e-1);
maxerr = max(error);
err = any(err);
data = [];

if opt.Display
    figure(9),clf
    subplot(122),hold on
    plot(r,P)
    plot(r,Result)
    Result = fitregmodel(S,K,r,L,'tikhonov',RegParam2,'Solver','fnnls');
    plot(r,Result)
    subplot(121),hold on
    plot(t,DipEvoFcn)
    plot(t,Filtered)
    plot(t,S)
end

end