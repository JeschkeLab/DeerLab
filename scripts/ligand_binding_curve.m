% Two-component titration
clc, clear, clf

% Chemistry Definition
%---------------------------------------
%Prepare equilibrium of type:
%   A + L <-> B
KD = 5.65;  % dissociation constant, uM
Ctot = 1; % total protein concentration, uM
L = [0.1 0.3 1 3 10 30 100 300]; % total ligand concentration, uM
nDataSets = numel(L);
% Calculate mole fractions a (for A = protein without ligand) and
% b (for B = protein with bound protein)
Kb = 1/KD;
for q = 1:nDataSets
    b_ = roots([Kb*Ctot -(Kb*L(q)+Kb*Ctot+1) Kb*L(q)]);
    b(q) = b_(b_<=1 & b_>=0);
end
a = 1-b;

% Preparation
%---------------------------------------
t = linspace(-0.1,3,200);
r = time2dist(t);
ra = 3; wa = 0.2;
rb = 3.5; wb= 0.3;
pmodel = [ra wa rb wb];
% Generate set of DEER traces for different ligand concentrations
for i = 1:nDataSets
    P = rd_twogaussian(r,[pmodel b(i)]);
    V{i} = dipolarsignal(t,r,P,'ModDepth',1,'noiselevel',0.05);
end
K = dipolarkernel(t,r);
[Ks{1:nDataSets}] = deal(K);


% Fitting
%---------------------------------------
p0 = [linspace(0,1,nDataSets) 2 0.5 3 0.5];
plow = [zeros(1,nDataSets) 1 0.1 1 0.1];
pup = [ones(1,nDataSets) 5 0.7 5 0.7];

model = @(r,p,idx)rd_twogaussian(r,[p(nDataSets+1:end) p(idx)]);
pfit = fitparamodel(V,model,r,Ks,p0,'Lower',plow,'Upper',pup);

%Get the fitted molar fractions
bfit = pfit(1:8);
afit = 1 - bfit;
%Control label switching
if mean(diff(bfit))<0
    tmp = afit;
    afit=bfit;
    bfit = tmp;
end

for i=1:nDataSets
    Pfit{i} = rd_twogaussian(r,[pfit(nDataSets+1:end) pfit(i)]);
    Vfit{i} = Ks{i}*Pfit{i};
end

%Fit the dissociation constant to the fitted molar fractions
Kbfit = fminsearch(@(Kb)norm(Kb*Ctot*bfit.^2 -(Kb*L+Kb*Ctot+1).*bfit + Kb*L)^2,1);
KDfit = 1/Kbfit;

% Plotting
%---------------------------------------
%Extrapolate the fitted molar fraction values for plotting
Lfine = logspace(min(log10(L)),max(log10(L)),50);
for q = 1:numel(Lfine)
    p = [Kbfit*Ctot -(Kbfit*Lfine(q)+Kbfit*Ctot+1) Kbfit*Lfine(q)];
    b_ = roots(p);
    b(q) = b_(b_<=1 & b_>=0);
end
a = 1-b;

%Plots
subplot(221)
plot(t,cell2mat(V)+(0:nDataSets-1),'k.',t,cell2mat(Vfit)+(0:nDataSets-1),'r','LineWidth',1.5)
axis tight,grid on
xlabel('time [\mus]'),ylabel('V_i(t)')
legend('Data','Fit')
subplot(222)
plot(t,cell2mat(Pfit)+(0:nDataSets-1),'k','LineWidth',1.5)
axis tight,grid on
xlabel('distance [nm]'),ylabel('P_i(r)')
legend('Fit')
subplot(2,2,[3 4])
plot(log(L),bfit,'b.',log(L),afit,'r.',log(Lfine),b,'b',log(Lfine),a,'r','LineWidth',1.5,'MarkerSize',25)
axis tight,grid on
xlabel('log(Ligand conc.)'),ylabel('Molar fraction')
title(sprintf('K_D = %.2f/%.2f',KDfit,KD))
legend('#1 Exp.','#2 Exp.','#1 Fit','#2 Fit')
