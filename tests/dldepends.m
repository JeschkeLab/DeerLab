function dadepends(FileMask)

DeerLabPath = fileparts(which('dadepends'));
DeerLabPath = DeerLabPath(1:end-numel('\tests'));
if isempty(DeerLabPath)
    error('DeerLab is not on the Matlab path!');
end

if nargin==0
    FileMask = '*';
else
    FileMask = [FileMask '*'];
end
fid = 1;
fprintf(fid,'=======================================================================\n');
fprintf(fid,'DeerLab Toolbox Dependency Analysis           %s\n(Matlab %s)\n',datestr(now),version);
fprintf(fid,'DeerLab location: %s\n',DeerLabPath);
fprintf(fid,'=======================================================================\n');

%Run through ./tests/
files = dir(fullfile(DeerLabPath,'tests',FileMask));
filenames = {files.name};
if ~isempty(filenames)
    if strcmp(filenames{1},'.')
        filenames = {filenames{3:end}};
    end
end
fprintf(fid,'-----------------------------------------------------------------------\n');
fprintf(fid,'Path: /tests/                                                          \n');
fprintf(fid,'-----------------------------------------------------------------------\n');
fprintf(fid,'Filename                              Dependencies                   \n');

for i=1:length(filenames)
    
    [~,output] = matlab.codetools.requiredFilesAndProducts({fullfile(DeerLabPath,'tests',filenames{i})});
    
    if ~isempty(output)
        output = {output.Name};
        string = output{1};
        for j=2:length(output)
            string = [string ', ' output{j}];
        end
        
        fprintf(fid,'%-36s  %s \n',filenames{i},string);
    end
end

%Run through ./functions/
files = dir(fullfile(DeerLabPath,'functions',FileMask));
filenames = {files.name};
if ~isempty(filenames)
    if strcmp(filenames{1},'.')
        filenames = {filenames{3:end}};
    end
end
fprintf(fid,'-----------------------------------------------------------------------\n');
fprintf(fid,'Path: /functions/                                                          \n');
fprintf(fid,'-----------------------------------------------------------------------\n');
fprintf(fid,'Filename                              Dependencies                   \n');

for i=1:length(filenames)
    
    [~,output] = matlab.codetools.requiredFilesAndProducts({fullfile(DeerLabPath,'/functions/',filenames{i})});
    
    if ~isempty(output)
        output = {output.Name};
        string = output{1};
        for j=2:length(output)
            string = [string ', ' output{j}];
        end
        
        fprintf(fid,'%-36s  %s \n',filenames{i},string);
    end
    
end


fprintf(fid,'=======================================================================\n');

end


