	function [ x,y,Z, VB ]= get_elexsys(file,stop_VB)
% function [ Z, { [ BP_table values] } ]= get_elexsys(file)
% Get Bruker Data Set
% file is the file name for the bruker data set
% Use the Bruker  file convention, and do not add extension names
% "file" . dsc and .dta are assumed
% Z contains the x and y values in a 2 column matrix
%  plot(Z(:,1),Z(:,2))  to  visualize results
% VB contains as ascii strings the Bruker Mnemonics followed by the values of the mnemonics
%  uses get_vb2(VB,'mnemonic') to get the value of a parameter by mnemonic
% based on routines by B. Robinson for CW EPR spectra in Bruker ESP format
% modified Mar 2002 for Bruker Elexsys / Linux files by DH and GJ
% and Jan 2004 for 1D and 2D pulse EPR spectra that may be complex by GJ 

if ~exist('stop_VB','var'),
    stop_VB = 100;
end;

x=[];
y=[];
Z =[];
VB=[];
N=1;
M=1;
%keyboard
fid = fopen( [ file '.DTA'], 'r','ieee-be.l64') ;% spc
if(fid == -1) return; end
[Z,N]= fread(fid,inf,'float64');
fclose(fid);
fid = fopen( [ file '.DSC'], 'r');% par
[S,sn]= fscanf(fid,'%c');
fclose(fid);

%Parse string S, the delimiter is ascii 10, this is the token delimiter
	D = setstr(10);
	VB=[];
    poi = 0;
	while(length(S) > 0 && poi < stop_VB)
        poi = poi + 1;
	[token,S] = strtok(S,D);
	VB = str2mat(VB,token);
	end
 %keyboard

 GST = get_vb2(VB,'XMIN'); GSI = get_vb2(VB,'XWID'); N = get_vb2(VB,'XPTS');
	dx = GSI/(N-1);
    x=GST + dx*[0:(N-1)];
    dim2=get_vb2(VB,'YTYP');
    if ~strncmp(dim2,'NODATA',3),
        ymin = get_vb2(VB,'YMIN'); 
        ywid = get_vb2(VB,'YWID'); 
        M = get_vb2(VB,'YPTS');
	    dy = ywid/(M-1);
        y=ymin + dy*[0:(M-1)];
    end;
    
%     disp(length(Z));
%     disp(N);
%     disp(M);
    
    if length(Z)==2*N*M,
        z0=Z;
        Z=zeros(N*M,1);
        for k=1:N*M,
            Z(k)=z0(2*k-1)+1i*z0(2*k);
        end;
    end;

    Z=reshape(Z,N,M);
    Z=permute(Z,[2,1]);
%     disp(size(Z));
