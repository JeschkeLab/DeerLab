function distr=Sphere_Surface(r0,par)
%
% Model library of DeerAnalysis2006: Sphere_Surface
%
% a distribution on a spherical surface for spheres with mean diameter
% ds and a Gaussian distribution of mean diameters with standard deviation s(ds) 
%
% (c) G. Jeschke, 2006
%
% PARAMETERS
% name    symbol default lower bound upper bound
% par(1)  ds     4.0     1.5         10             mean radius of the sphere
% par(2)  s(ds)  0.3     0.05        5              standard dev. of sphere radius

rstart=min(r0);
gauss2=(r0-par(1)*ones(size(r0)))/(sqrt(2)*par(2));
gauss2=exp(-gauss2.^2);
minr=par(1)-3*par(2);
minarg=r0-minr*ones(size(r0));
maxr=par(1)+3*par(2);
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
distr=triangle;



