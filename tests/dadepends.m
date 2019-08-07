function dadepends(FileMask)

DeerAnalysisPath = fileparts(which('DeerAnalysis'));
if isempty(DeerAnalysisPath)
    error('DeerAnalysis is not on the Matlab path!');
end

if nargin==0
    FileMask = '*';
else
    FileMask = [FileMask '*'];
end
fid = 1;
fprintf(fid,'=======================================================================\n');
fprintf(fid,'DeerAnalysis Toolbox Dependency Analysis           %s\n(Matlab %s)\n',datestr(now),version);
fprintf(fid,'DeerAnalysis location: %s\n',DeerAnalysisPath);
fprintf(fid,'=======================================================================\n');

%Run through ./tests/
files = dir(fullfile(DeerAnalysisPath,'tests',FileMask));
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
    
    [~,output] = matlab.codetools.requiredFilesAndProducts({fullfile(DeerAnalysisPath,'tests',filenames{i})});
    
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
files = dir(fullfile(DeerAnalysisPath,'functions',FileMask));
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
    
    [~,output] = matlab.codetools.requiredFilesAndProducts({fullfile(DeerAnalysisPath,'/functions/',filenames{i})});
    
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


