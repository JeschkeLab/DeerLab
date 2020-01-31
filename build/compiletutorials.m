
clc

files = dir('../tutorials/**/*.mlx');
sourcenames = {files.name};
pdfnames = cellfun(@(S)[S(1:end-3) 'pdf'],sourcenames,'UniformOutput',false);
htmlnames = cellfun(@(S)[S(1:end-3) 'html'],sourcenames,'UniformOutput',false);
mnames = cellfun(@(S)[S(1:end-3) 'm'],sourcenames,'UniformOutput',false);

paths = {files.folder};
sourcefiles = cellfun(@(x,y)fullfile(x,y),paths,sourcenames,'UniformOutput',false);
pdffiles = cellfun(@(x,y)fullfile(x,y),paths,pdfnames,'UniformOutput',false);
htmlfiles = cellfun(@(x,y)fullfile(x,y),paths,htmlnames,'UniformOutput',false);
mfiles = cellfun(@(x,y)fullfile(x,y),paths,mnames,'UniformOutput',false);

%Convert live script into other formats
for i=1:numel(pdfnames)
    fprintf('Conversion: mlx -> pdf (file %i of %i)...',i,numel(sourcefiles))
    matlab.internal.liveeditor.openAndConvert(sourcefiles{i},pdffiles{i});
    fprintf(' complete \n')
end

for i=1:numel(htmlnames)
    fprintf('Conversion: mlx -> html (file %i of %i)...',i,numel(sourcefiles))
    matlab.internal.liveeditor.openAndConvert(sourcefiles{i},htmlfiles{i});
    fprintf(' complete \n')
end

for i=1:numel(mnames)
    fprintf('Conversion: mlx -> m (file %i of %i)...',i,numel(sourcefiles))
    matlab.internal.liveeditor.openAndConvert(sourcefiles{i},mfiles{i});
    fprintf(' complete \n')
end


if ispc
    
    %Synchronize the AWS bucket
    system('python synchS3.py --keyfile %USERPROFILE%\.ssh\aws_access_keys.txt --directory "../tutorials" --bucket deertutorials');
    
    %Delete files
    for i=1:numel(pdfnames)
        fprintf('Deleting: pdf (file %i of %i)...',i,numel(sourcefiles))
        delete(pdffiles{i});
        fprintf(' complete \n')
    end
    
    for i=1:numel(htmlnames)
        fprintf('Deleting: html (file %i of %i)...',i,numel(sourcefiles))
        delete(htmlfiles{i});
        fprintf(' complete \n')
    end
    
    for i=1:numel(htmlnames)
        fprintf('Deleting: m (file %i of %i)...',i,numel(sourcefiles))
        delete(mfiles{i});
        fprintf(' complete \n')
    end
    
end
