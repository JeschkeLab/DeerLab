function [err,data,maxerr] = test(opt,olddata)


    t = linspace(-1,5,200);
    r= time2dist(t);
    P = rd_onegaussian(r,[5,0.2]);
    B = td_exp(t,0.3);
    TrueOffset = 1e8;
    noiselevel = 0.05;
    V = dipolarsignal(t,r,P,'Moddepth',0.25,'Background',B,...
        'Offset',TrueOffset,'noiselevel',noiselevel);


    [Vc,Offset] = correctscale(V,t);
    
    
    err  = abs(Offset - TrueOffset)/TrueOffset>noiselevel;
    data = [];
    maxerr = max( abs(Offset - TrueOffset)/TrueOffset);
    
    
    if opt.Display
       figure(8)
       plot(t,Vc)
        
    end



end