function [r,pcf,pop]=conformers2PCF(fname,Ec,Dc,T)
% 
% function [r,pcf]=conformers2PCF(fname,Ec,Dc)
%
% Computes a pair correlation function (distance distribution) from a list
% of conformers given in the ASCII file fname (extension .dat is assumed)
% populations of the conformers according to a Boltzmann distribution are 
% also computed
%
% Ec    columns in the ASCII file that contain contributions to the energy 
%       of the conformer, needs to be in kJ/mol (can be a single number or
%       vector)
% Dc    columns in the ASCI file that contain distances that have to be
%       averaged for each conformer (can be a single number or vector)
% T     temperature (optional), default value is 298 K
%
% r     distance axis
% pcf   pair correlation function, normalized to unity sum
% pop   populations of the conformers, normalized to unity sum
%
% Example:
% an ASCII file 'testfile.dat' for a nitroxide biradical may start as follows
% (note the percent (comment) character at the beginning of the table headline 
% which is required by MATLAB):
% %E (kJ/mol)    r(N-N) (A)      r(O-O) (A)
% 791.642914	 31.4779270	     31.2486643	
% 797.332140	 40.4554192	     41.5977748
% ...
%
% the proper call for the distance distribution at 244 K (glass transition
% temperature of o-terphenyl) would be:
% [r,pcf,pop]=conformers2PCF('testfile',1,[2 3],244);
%
% (c) G. Jeschke, 2007

R=8.314; % universal gas constant in J mol^{-1} K^{-1}

% if no temperature is given, set to room temperature default
if nargin<4,
    T=298;
end;

data=load([fname '.dat']);
[m,n]=size(data);
E=sum(data(:,Ec),2)/length(Ec); % energy average
E=E-min(E)*ones(size(E)); % put energy minimum to zero
D=sum(data(:,Dc),2)/length(Dc); % distance average

dr=max(D)-min(D); % width of distance range
r=linspace(min(D)-0.2*dr,max(D)+0.2*dr,round(m/5)); % distance axis
ra=min(r);
re=max(r);
rn=length(r);
ddr=r(2)-r(1);
pcf=zeros(1,round(m/5)); % initialize PCF
pop=zeros(1,m); % initialize populations

for k=1:m,
    rpoi=1+floor((rn-1)*(D(k)-ra)/(re-ra));
    delr=D(k)-r(rpoi);
    pop(k)=exp(-1.0e3*E(k)/(R*T));
    pcf(rpoi)=pcf(rpoi)+pop(k)*(1-delr/ddr);
    pcf(rpoi+1)=pcf(rpoi+1)+pop(k)*delr/ddr;
end;
z=sum(pop); % partition function
pop=pop/z;
pcf=pcf/z;
