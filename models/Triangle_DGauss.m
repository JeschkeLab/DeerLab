function [deer,distr]=Triangle_DGauss(r0,t0,par)
%
% Model library of DeerAnalysis2008: triangle with Gaussian distribution of
% vertex positions
%
% two Gaussian peaks with mean distance <r> and width (standard
% deviation) s(r)
% (c) G. Jeschke, 2009
%
% #extended# denotes a model that provides both distribution and deer trace
% #enable# 5 only the five first parameters are fitted by default
% 
% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  <rv>   2.5     0.5         10         mean distance from C3 axis
% par(2)  s(v)   0.5     0.02        5          std. dev. of vertex position
% par(3)  p1     0.5     0.05        0.95       population of 1st component
% par(4)  <rv2>  1.5     0.5         10         2nd mean distance
% par(5)  s(v2)  0.5     0.02        5          std. dev of 2nd vertex pos.
% par(6)  Delta    1     0.1         1          total modulation depth
% par(7)  nmc    2500    1000        100000     number of Monte Carlo trials



rv=par(1);
sv=par(2);
pop1=par(3);
pop2=1-pop1;
rv2=par(4);
sv2=par(5);
Delta=par(6);
lambda=1-sqrt(1-Delta);
nmc=par(7);
distr=zeros(size(r0));
pair=zeros(size(t0));
triple=zeros(size(t0));
ra=r0(1);
re=r0(end);
nr=length(r0)-1;

for k=1:nmc
    % first component
    fi1=2*pi*rand; % uniform distribution of polar angle of vertex shift vector
    dv=sv*randn; % normal distribution of vertex distance from mean position
    p1=[0,rv]+[cos(fi1)*dv,sin(fi1)*dv]; % vertex 1 coordinates
    fi1=2*pi*rand; % uniform distribution of polar angle of vertex shift vector
    dv=sv*randn; % normal distribution of vertex distance from mean position
    p2=[sqrt(3)*rv/2,-rv/2]+[cos(fi1)*dv,sin(fi1)*dv]; % vertex 2 coordinates
    fi1=2*pi*rand; % uniform distribution of polar angle of vertex shift vector
    dv=sv*randn; % normal distribution of vertex distance from mean position
    p3=[-sqrt(3)*rv/2,-rv/2]+[cos(fi1)*dv,sin(fi1)*dv]; % vertex 3 coordinates
    d12=norm(p1-p2); % side length 1,2
    poi12=1+round(nr*(d12-ra)/(re-ra)); % pointer into distance distribution
    if poi12>0 && poi12<=nr+1
        distr(poi12)=distr(poi12)+pop1;
    end
    d13=norm(p1-p3); % side length 1,3
    poi13=1+round(nr*(d13-ra)/(re-ra)); % pointer into distance distribution
    if poi13>0 && poi13<=nr+1
      distr(poi13)=distr(poi13)+pop1;
    end
    d23=norm(p2-p3); % side length 1,2
    poi23=1+round(nr*(d23-ra)/(re-ra)); % pointer into distance distribution
    if poi23>0 && poi23<=nr+1
      distr(poi23)=distr(poi23)+pop1;
    end
    fi2=2*pi*rand;
    cth2=rand;
    sth2=sqrt(1-cth2^2);
    B0=[cos(fi2)*sth2,sin(fi2)*sth2];
    cos1=sum(((p1-p2)/d12).*B0);
    wd1=2*pi*(1-3*cos1^2)*52.04/d12^3;
    mod1=cos(wd1*t0);
    cos2=sum(((p1-p3)/d13).*B0);
    wd2=2*pi*(1-3*cos2^2)*52.04/d13^3;
    mod2=cos(wd2*t0);
    cos3=sum(((p2-p3)/d23).*B0);
    wd3=2*pi*(1-3*cos3^2)*52.04/d23^3;
    mod3=cos(wd3*t0);
    pair=pair+pop1*(mod1+mod2+mod3)/3;
    triple=triple+pop1*(mod1.*mod2+mod1.*mod3+mod2.*mod3)/3;
    % second component
    fi1=2*pi*rand; % uniform distribution of polar angle of vertex shift vector
    dv=sv2*randn; % normal distribution of vertex distance from mean position
    p1=[0,rv2]+[cos(fi1)*dv,sin(fi1)*dv]; % vertex 1 coordinates
    fi1=2*pi*rand; % uniform distribution of polar angle of vertex shift vector
    dv=sv2*randn; % normal distribution of vertex distance from mean position
    p2=[sqrt(3)*rv2/2,-rv2/2]+[cos(fi1)*dv,sin(fi1)*dv]; % vertex 2 coordinates
    fi1=2*pi*rand; % uniform distribution of polar angle of vertex shift vector
    dv=sv2*randn; % normal distribution of vertex distance from mean position
    p3=[-sqrt(3)*rv2/2,-rv2/2]+[cos(fi1)*dv,sin(fi1)*dv]; % vertex 3 coordinates
    d12=norm(p1-p2); % side length 1,2
    poi12=1+round(nr*(d12-ra)/(re-ra)); % pointer into distance distribution
    if poi12>0 && poi12<=nr+1
        distr(poi12)=distr(poi12)+pop2;
    end
    d13=norm(p1-p3); % side length 1,3
    poi13=1+round(nr*(d13-ra)/(re-ra)); % pointer into distance distribution
    if poi13>0 && poi13<=nr+1
      distr(poi13)=distr(poi13)+pop2;
    end
    d23=norm(p2-p3); % side length 1,2
    poi23=1+round(nr*(d23-ra)/(re-ra)); % pointer into distance distribution
    if poi23>0 && poi23<=nr+1
      distr(poi23)=distr(poi23)+pop2;
    end
    fi2=2*pi*rand;
    cth2=rand;
    sth2=sqrt(1-cth2^2);
    B0=[cos(fi2)*sth2,sin(fi2)*sth2];
    cos1=sum(((p1-p2)/d12).*B0);
    wd1=2*pi*(1-3*cos1^2)*52.04/d12^3;
    mod1=cos(wd1*t0);
    cos2=sum(((p1-p3)/d13).*B0);
    wd2=2*pi*(1-3*cos2^2)*52.04/d13^3;
    mod2=cos(wd2*t0);
    cos3=sum(((p2-p3)/d23).*B0);
    wd3=2*pi*(1-3*cos3^2)*52.04/d23^3;
    mod3=cos(wd3*t0);
    pair=pair+pop2*(mod1+mod2+mod3)/3;
    triple=triple+pop2*(mod1.*mod2+mod1.*mod3+mod2.*mod3)/3;
end

distr=distr/sum(distr); % normalization
deer=2*lambda*(1-lambda)*pair+lambda^2*triple;
deer=deer/max(deer);

