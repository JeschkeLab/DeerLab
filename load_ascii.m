function dset=load_ascii(inname,header_lines)
%
% Load data array from ASCII file with known number of header lines
%
% inname        filename (pathname) of ASCII data file
% header_lines  number of header lines before data section
%
% dset          output data array
%
% limitation: 2048 data rows, 10 columns
%
% G. Jeschke, 2002
%
% for comments on implementation, see sister subroutine in DeerTrafo
%

rfile=fopen(inname,'r');

nl=0;
while nl<header_lines
   fgetl(rfile);
   nl=nl+1;
end
dset=zeros(2048,10);
line=' ';
poi=0;
if nl==0, line=' '; end

while ~isempty(line)
   line=fgetl(rfile);
   if ~isstr(line) || length(line)<2, break,
   else
       vec=str2num(line);
       if poi==0, columns=length(vec); end
       if length(vec)==columns
           poi=poi+1;
           dset(poi,1:length(vec))=vec;
       end
   end
end

dset=dset(1:poi,1:columns);

fclose(rfile);