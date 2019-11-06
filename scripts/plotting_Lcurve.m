%======================================================================
% DeerAnalyis2
% Example: Selection methods and the L-curve
% Computing different alpha-selection methods and plottig them on the 
% L-curve for comparison.
%=======================================================================

clear,clc,clf

%Prepare signal and kernel
t = linspace(-0.1,6,200);
r = linspace(2,8,200);
V = (dipolarsignal(t,r,rd_onegaussian(r,[4 0.2]),'noiselevel',0.1));
K = dipolarkernel(t,r);

%Get the optimal regularization parameters
[alpha,~,alphas,Res,Pen] = selregparam(V,K,r,'tikh',{'aic','mcl','gml','lr'});

%Find the positions of the different optimal regularizfation parameters 
%on the L-curve
AICidx = alphas==alpha(1);
BICidx = alphas==alpha(2);
GMLidx = alphas==alpha(3);
LRidx  = alphas==alpha(4);

%Construct the vector of the L-curve
logRes = log(Res.^2);
logPen = log(Pen.^2);

%Plot the L-cruve and the optimal parameters
clf,hold on
plot(logRes,logPen,'LineWidth',1.5)
plot(logRes(AICidx),logPen(AICidx),'.','MarkerSize',25)
plot(logRes(BICidx),logPen(BICidx),'.','MarkerSize',25)
plot(logRes(GMLidx),logPen(GMLidx),'.','MarkerSize',25)
plot(logRes(LRidx),logPen(LRidx),'.','MarkerSize',25)
legend('L-curve','AIC','BIC','GML','LR')
grid on,axis tight, box on
xlabel('log(||KP - V)||^2')
ylabel('log(||LP)||^2')
set(gca,'FontSize',13)

