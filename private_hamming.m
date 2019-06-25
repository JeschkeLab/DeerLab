function spcn = private_hamming(spc,r)
%HAMMING	Hamming apodization of a spectrum, renamed private_hamming 
%           to avoid conflicts with the Signal Processing Toolbox 
%	spcn = private_hamming(spc,r)
%
%  r part at which hamming window has decayed to zero
%	

%	Copyright (c) 1998 by Gunnar Jeschke

[m,n]=size(spc);
spcn=spc;
rad=round(r*n);
arg=linspace(0,pi,rad);
hamm=0.54*ones(1,rad)+0.46*cos(arg);
if rad>n, hamm=hamm(1,1:n); end;
if rad<n, hamm=[hamm zeros(1,n-rad)]; end;
for k=1:m,
   spcn(k,:)=hamm.*spcn(k,:);
end;

		

