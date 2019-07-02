function distr=WLC_rigid(r0,par)
%
% Model library of DeerAnalysis2006: WLC_rigid
%
% Worm-like chain model near the rigid limit
% J. Wilhelm, E. Frey, Phys. Rev. Lett. 77(12), 2581-2584 (1996)
% equations (3) and (4)
% L is the length of the worm-like chain
% Lp is the persistence length
%
% (c) G. Jeschke, 2006
%
% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  L      3.7     1.5         10             length of the worm-like chain
% par(2)  Lp     10      2           100            persistence length

L=par(1);
Lp=par(2);
kappa=Lp/L;

rn=r0/L;
distr=zeros(size(r0));
terms=2;
for pp=1:length(rn)
    G=0;
    crit=kappa*(1-rn(pp));
    if crit>0.2
        fac=2*kappa/(4*pi);
        for k=1:terms
            G=G+fac*pi^2*k^2*(-1)^(k+1)*exp(-kappa*pi^2*k^2*(1-rn(pp)));
        end
    elseif crit>0
        fac=kappa/(4*pi*2*sqrt(pi));
        for l=1:terms
            harg=(l-1/2)/sqrt(kappa*(1-rn(pp)));
            h2=4*harg^2-2;
            G=G+fac*1/(kappa*(1-rn(pp)))^(3/2)*exp(-(l-1/2)^2/(kappa*(1-rn(pp))))*h2;
        end
    end
    distr(pp)=G;
end

distr = distr/sum(distr);

