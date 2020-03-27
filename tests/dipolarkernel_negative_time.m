function [pass,maxerr] = test(opt)

% Check that kernel is constructed properly for negative times using fresnel method

tneg = linspace(-5,0,50);
tpos = linspace(0,5,50);
r = 3;
Kneg = dipolarkernel(tneg,r);
Kpos = dipolarkernel(tpos,r);

delta = abs(Kneg - flipud(Kpos));

% Pass: the negative and positive kernel parts are identical
pass = all(delta(:) < 1e-12);

maxerr = max(delta(:));

%Plot if requested
if opt.Display
    plot(tneg,Kneg,tpos,Kpos);
    legend('negative','positive');
    xlabel('t [\mus]')
    ylabel('V(t)')
    grid on, axis tight, box on
end

end