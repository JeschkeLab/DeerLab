	function [ value , vstr] = get_vb2(VB,mnem)
%  [ value , vstring ] = get_vb2(VB,'mnem')
% VB contains the Bruker menmonic list with the values from a data set
% menm is a known bruker menmonic. e.g. GSI is starting field in gauss
% also get the value as an uncoverted string if desired
% if vstr cannot be converted to a numeric value
% then value and vstr are the same, except that all spaces are removed from value,
% and are both strings
% this is useful for comments etc.
% value and vstr come back blank if not found
	[m,n] = size(VB);
	value = [];
	vstr =[];
	k =1;
	while( (k<m) & ~ strcmp( strtok(VB(k,:)),mnem ) )
	k=k+1;, end
	if (k<m) [t,vstr] = strtok(VB(k,:));
		value = str2num(vstr); end
	if(isempty(value)),
        vstr0=vstr;
        vstr=[];
        for k=1:length(vstr0), if ~isspace(vstr0(k)), vstr=[vstr vstr0(k)]; end; end;
        value = vstr;, 
    end

	




	