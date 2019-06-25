function netsets = get_netset_definitions(fname)

fid=fopen(fname);
if fid==-1
    netsets(1).name = '';
    netsets(1).directory = '';
    netsets(1).version = [];
    netsets(1).url = '';
    return;
end;

nl = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~isempty(tline)
        k = strfind(tline,'%'); % remove comments
        if ~isempty(k)
            if k(1)>1
                tline = tline(1:k(1)-1);
            else
                tline = '%';
            end;
        end;
        myline = textscan(tline,'%s','Delimiter','|');
        args=myline{1};
        nl = nl + 1;
        netsets(nl).name = char(args(1));
        netsets(nl).directory = char(args(2));
        netsets(nl).version = str2double(char(args(3)));
        netsets(nl).url = char(args(4));        
    end
end

fclose(fid);