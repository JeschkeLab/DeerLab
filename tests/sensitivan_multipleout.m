function [pass,maxerr] = test(opt)

M = 200;
t = linspace(0,4,M);
r = time2dist(t);
B = td_exp(t,0.3);
P = rd_onegaussian(r,[4,0.3]);
V = dipolarsignal(t,r,P,'noiselevel',0.05,'ModDepth',0.3,'Background',B);

Parameters.regparam = linspace(10,50,2);
Parameters.validationnoise = linspace(0.01,0.1,2);

if opt.Display
    f = figure(1); clf;AxisHandle = axes(f);
else
    AxisHandle = [];
end

fcnHandle = @(param)myfitting(param,t,r,V);

stats  = sensitivan(fcnHandle,Parameters,'AxisHandle',AxisHandle);

err(1) = ~isstruct(stats);
err(2) = length(stats)~=2;
pass = all(err);
 
maxerr = 0;

if opt.Display
    cla
    subplot(121)
    hold on
    plot(r,meansOut{1},'b','LineWidth',1)
    f = fill([r fliplr(r)] ,[Upper{1} fliplr(Lower{1})],...
        'b','LineStyle','none');
    f.FaceAlpha = 0.5;
    subplot(122)
    axis tight, grid on, box on
    hold on
    plot(t,meansOut{2},'b','LineWidth',1)
    f = fill([t fliplr(t)] ,[Upper{2}  fliplr(Lower{2})],...
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