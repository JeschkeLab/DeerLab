	function [ x,y,z, VB ]= Get_Brux(file)
% function [ Z, { [ BP_table values] } ]= Get_Brux(file)
% Get Bruker Data Set
% file is the file name for the bruker data set
% Use the Bruker  file convention, and do not add extension names
% "file" . par and .spc are assumed
% Z contains the x and y values in a 2 column matrix
%  plot(Z(:,1),Z(:,2))  to  visualize results
% VB contains as ascii strings the Bruker Mnemonics followed by the values of the mnemonics
%  use get_vb(VB,'mnemonic') to get the value of a parameter by mnemonic

z =[];
%keyboard
fid = fopen( [ file '.spc'], 'r') ;
if(fid == -1) return;, end
[z,N]= fread(fid,inf,'float');
fclose(fid);
fid = fopen( [ file '.par'], 'r');
[S,sn]= fscanf(fid,'%c');
fclose(fid);

%Parse string S, the delimeter is ascii 13, this is the token delimiter
	D = setstr(13);
	VB=[];
	while(length(S) > 0 )
	[token,S] = strtok(S,D);
	VB = str2mat(VB,token);
	end
 %keyboard


 	GST = get_vb2(VB,'GST'), GSI = get_vb2(VB,'GSI'), N = get_vb2(VB,'XPLS'),
	dx = GSI/N;
    x=GST+dx*[0:(N-1)];
    y=[];
%	keyboard