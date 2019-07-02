function distr=Chechik_2(r0,par)
%
% Model library of DeerAnalysis2006: Chechik_2
%
% a Gaussian peak with mean distances <r1> and standard deviation s(r1)
% plus a distribution on a spherical surface for spheres with mean diameter
% ds and a Gaussian distribution of mean diameters with standard deviation s(ds) 
% the peaks have relative integral intensities
% p1 and p2=1-p1
%
% (c) G. Jeschke, 2006
%
% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  <r1>   2.0     1.0         10
% par(2)  s(r1)  0.3     0.02        5
% par(3)  p1     0.25    0           1
% par(4)  ds     4.0     1.5         10
% par(5)  s(ds)  0.3     0.02        5

rstart=min(r0);
gauss1=(r0-par(1)*ones(size(r0)))/(sqrt(2)*par(2));
gauss1=exp(-gauss1.^2);
intg1=sum(gauss1);
gauss2=(r0-par(4)*ones(size(r0)))/(sqrt(2)*par(5));
gauss2=exp(-gauss2.^2);
minr=par(4)-3*par(5);
minarg=r0-minr*ones(size(r0));
maxr=par(4)+3*par(5);
maxarg=r0-maxr*ones(size(r0));
[~,mapoi]=min(abs(maxarg));
[~,mipoi]=min(abs(minarg));
triangle=zeros(size(r0));
for k=mipoi:mapoi
    rac=r0(k);
    he=1/rac; % scale amplitude so that area is constant
    he=he*rac^2; % rescale amplitude, considering that number of labels scales with sphere surface
    ha=rstart*he/rac; % amplitude at first point
    triarg=r0-r0(k)*ones(size(r0));
    [~,poi]=min(abs(triarg));
    triarg=linspace(ha,he,poi);
    triangle(1:poi)=triangle(1:poi)+gauss2(k)*triarg;
end
triangle=triangle*intg1/sum(triangle);
distr=par(3)*gauss1+(1-par(3))*triangle;
distr = distr/mean(distr);
