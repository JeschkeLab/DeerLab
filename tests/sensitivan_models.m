function [err,data,maxerr] = test(opt,olddata)

M = 200;
t = linspace(0,4,M);
r = time2dist(t);
B = td_exp(t,0.3);
P = rd_onegaussian(r,[4,0.3]);
V = dipolarsignal(t,r,P,'noiselevel',0.02,'ModDepth',0.3,'Background',B);

Parameters.regparam = linspace(10,50,3);
Parameters.Bmodel = {@td_exp,@td_strexp};

if opt.Display
    f = figure(1); clf;AxisHandle = axes(f);
else
    AxisHandle = [];
end

fcnHandle = @(param)myfitting(param,t,r,V);

[meansOut,Upper,Lower] = sensitivan(fcnHandle,Parameters,'AxisHandle',AxisHandle);

err(1) = ~iscell(meansOut);
err(2) = length(meansOut)~=2;
err = any(err);
data = [];
maxerr = 0;

if opt.Display
    cla
    subplot(121)
    hold on
    plot(r,meansOut{1},'b','LineWidth',1)
    f = fill([r fliplr(r)] ,[Upper{1} fliplr(Lower{1})],...
        'b','LineStyle','none');
    f.FaceAlpha = 0.5;
    axis tight, grid on, box on
    subplot(122)
    hold on
    plot(t,meansOut{2},'b','LineWidth',1)
    f = fill([t fliplr(t)] ,[Upper{2} fliplr(Lower{2})],...
        'b','LineStyle','none');
    f.FaceAlpha = 0.5;
    axis tight, grid on, box on
end

end


function [Pfit,Bfit] = myfitting(param,t,r,V)

[Bfit,lambdafit] = fitbackground(V,t,param.Bmodel);

K = dipolarkernel(t,r,lambdafit,Bfit);
regparam = param.regparam;
Pfit = fitregmodel(V,K,r,'tikh',regparam);

end