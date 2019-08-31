clc
%Load experimental data
[Expt,ExpS] = eprload('CT_DEER_mix_28_36.DSC');

%Correct phase and zero-time
S = correctPhase(ExpS);
[S,t] = correctZeroTime(S,Expt);
%Convert to us and truncate
t = t/1000;
[~,ZeroTimePos] = min(abs(t));
S = S(ZeroTimePos:end);
t = t(ZeroTimePos:end);

%Fit background
FitStartPos = 20;
Data2fit = S(FitStartPos:end);
Fitt = t(FitStartPos:end);
B = fitB(Data2fit,t,Fitt,'exponential');

%Ghost distance suppression
S = supressGhostDistances(S,2);

%Perform partial background correction
S = S./sqrt(B);

%Prepare regularization
r = time2dist(t);
K = getK(t,r,B);
RegMatrix = getRegMatrix(length(K),2);
RegParamRange = getRegParamRange(K,RegMatrix);
RegParam = selectRegParam(RegParamRange,S,K,RegMatrix,{'gml','lr'});

%Run regularization
Distribution1 = regularize(S,K,RegMatrix,'tikhonov',RegParam(1),'Solver','fnnls');
Distribution2 = regularize(S,K,RegMatrix,'tikhonov',RegParam(2),'Solver','fnnls');

figure(1),clf

subplot(121),hold on
plot(t,S,'k')
plot(t,sqrt(B))
plot(t,K*Distribution1,'r')
plot(t,K*Distribution2,'b')
xlabel('Time [\mus]')
ylabel('Intensity [a.u.]')
legend('Exp','Bckg','GML','L-Curve')

subplot(122),hold on
plot(r,Distribution1,'r')
plot(r,Distribution2,'b')
xlabel('Distance [nm]')
ylabel('P(r)')
legend('GML','L-Curve')
