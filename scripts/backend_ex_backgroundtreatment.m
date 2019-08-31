clc
%Load experimental data
[Expt,ExpS] = eprload('CT_DEER_mix_28_36.DSC');
Noise = rand(400,1);
Noise = Noise - mean(Noise);
Noise = 0.01*Noise/max(Noise)';
ExpS = ExpS/max(ExpS);
ExpS = ExpS + Noise;

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
B = fitB(Data2fit,t,Fitt);

%Ghost distance suppression
S = supressGhostDistances(S);

%Perform partial background correction
S2Reg = S./sqrt(B);

%Prepare regularization
r = time2dist(t);
K = getK(t,r,B);
RegMatrix = getRegMatrix(length(K),2);
RegParamRange = getRegParamRange(K,RegMatrix);
RegParam = selectRegParam(RegParamRange,S2Reg,K,RegMatrix,{'gml'});

%Run regularization
Distribution1 = regularize(S2Reg,K,RegMatrix,'tikhonov',RegParam(1),'Solver','fnnls');

S2Reg2 = S;

%Prepare regularization
K = getK(t,r,B,'KBType','full');
RegParamRange = getRegParamRange(K,RegMatrix);
RegParam = selectRegParam(RegParamRange,S2Reg2,K,RegMatrix,{'gml'});

Distribution2 = regularize(S2Reg2,K,RegMatrix,'tikhonov',RegParam(1),'Solver','fnnls');

ModDepth = 1/B(1) - 1;
S2Reg3 = S./B;
S2Reg3 = (S2Reg3 - (1-ModDepth))/ModDepth - 1;
% S2Reg3 = S2Reg3/S2Reg3(1);
%Prepare regularization
K = getK(t,r);
RegParamRange = getRegParamRange(K,RegMatrix);
RegParam = selectRegParam(RegParamRange,S2Reg3,K,RegMatrix,{'gml'});

Distribution3 = regularize(S2Reg3,K,RegMatrix,'tikhonov',RegParam(1),'Solver','fnnls');


figure(1),clf

subplot(121),hold on
plot(t,S2Reg3,'k')
plot(t,K*Distribution1)
plot(t,K*Distribution2)
plot(t,K*Distribution3)

xlabel('Time [\mus]')
ylabel('Intensity [a.u.]')
legend('Sqrt(B)','B','division')

subplot(122),hold on
plot(r,Distribution1)
plot(r,Distribution2)
plot(r,Distribution3)

xlabel('Distance [nm]')
ylabel('P(r)')
legend('Sqrt(B)','B','division')
