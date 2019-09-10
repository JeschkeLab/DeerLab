function [err,data,maxerr] = test(opt,olddata)

M = 200;
t = linspace(0,4,M);
r = time2dist(t);
B = td_exp(t,0.3);
P = rd_onegaussian(r,[4,0.3]);
V = dipolarsignal(t,r,P,'noiselevel',0.05,'ModDepth',0.3,'Background',B);

Parameters.regparam = linspace(10,50,5);
Parameters.validationnoise = linspace(0.01,0.1,2);

if opt.Display
    f = figure(1); clf;AxisHandle = axes(f);
else
    AxisHandle = [];
end

fcnHandle = @(param)myfitting(param,t,r,V);

[meansOut,stdsOut] = validate(fcnHandle,Parameters,'AxisHandle',AxisHandle);

err(1) = ~iscell(meansOut);
err(2) = length(meansOut)~=2;
err = any(err);
data = [];
maxerr = 0;

if opt.Display
    cla
    subplot(121)
    meanOut = meansOut{1};
    stdOut = stdsOut{1};
    hold on
    plot(r,meanOut,'b','LineWidth',1)
    f = fill([r fliplr(r)] ,[meanOut+stdOut fliplr(meanOut-stdOut)],...
        'b','LineStyle','none');
    f.FaceAlpha = 0.5;
    axis tight, grid on, box on
    subplot(122)
    meanOut = meansOut{2};
    stdOut = stdsOut{2};
    hold on
    plot(t,meanOut,'b','LineWidth',1)
    f = fill([t fliplr(t)] ,[meanOut+stdOut fliplr(meanOut-stdOut)],...
        'b','LineStyle','none');
    f.FaceAlpha = 0.5;
    axis tight, grid on, box on
end

end


function [Pfit,Bfit] = myfitting(param,t,r,V)

V = V + whitegaussnoise(length(V), param.validationnoise);

[Bfit,lambdafit] = fitbackground(V,t,@td_exp);

order = 2;
K = dipolarkernel(t,r,lambdafit,Bfit);
regparam = param.regparam;
Pfit = fitregmodel(V,K,r,'tikh',regparam);

end