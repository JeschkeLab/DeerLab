function [err,data,maxerr] = test(opt,olddata)

%My script
M = 200;
t = linspace(0,4,M);
r = time2dist(t);

P = rd_onegaussian(r,[4,0.3]);
K = dipolarkernel(t,r);
S = K*P;
validationnoise = 0;
S = dipolarsignal(t,r,P,'noiselevel',0.05);

S = S + whitegaussnoise(M,validationnoise);

order = 2;
L = regoperator(M,order);
regparam = 5;
Pfit = fitregmodel(S,K,r,L,'tikh',regparam);

Parameters(1).name = 'regparam';
Parameters(1).values = linspace(0.1,1,5);

Parameters(2).name = 'validationnoise';
Parameters(2).values = linspace(0.01,0.1,2);

if opt.Display
    f = figure(1); clf;AxisHandle = axes(f);
else
    AxisHandle = [];
end

[meanOut,stdOut] = validate('Pfit',Parameters,'AxisHandle',AxisHandle);

err = any(abs(P - meanOut)>5e-1);
data = [];
maxerr = max(abs(P - meanOut));

if opt.Display
        cla
        hold on
        plot(r,P,'k')
        plot(r,Pfit,'r--')
        plot(r,meanOut,'b','LineWidth',1)
        f = fill([r fliplr(r)] ,[meanOut.'+stdOut.' fliplr(meanOut.'-stdOut.')],...
            'b','LineStyle','none');
        f.FaceAlpha = 0.5;
        hold('off')
        axis('tight')
        grid('on')
        box('on')
end

end