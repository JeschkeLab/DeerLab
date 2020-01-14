

files = dir('../tutorials/**/*.mlx');
sourcenames = {files.name};
pdfnames = cellfun(@(S)[S(1:end-3) 'pdf'],sourcenames,'UniformOutput',false);
htmlnames = cellfun(@(S)[S(1:end-3) 'html'],sourcenames,'UniformOutput',false);
mnames = cellfun(@(S)[S(1:end-3) 'm'],sourcenames,'UniformOutput',false);

paths = {files.folder};
sourcefiles = cellfun(@(x,y)fullfile(x,y),paths,sourcenames,'UniformOutput',false);
pdffiles = cellfun(@(x,y)fullfile(x,y),paths,pdfnames,'UniformOutput',false);
htmlfiles = cellfun(@(x,y)fullfile(x,y),paths,pdfnames,'UniformOutput',false);
mfiles = cellfun(@(x,y)fullfile(x,y),paths,pdfnames,'UniformOutput',false);


for i=1:numel(pdfnames)
    fprintf('Conversion: mlx -> pdf (file %i of %i)...',i,numel(sourcefiles))
    matlab.internal.liveeditor.openAndConvert(sourcefiles{i},pdffiles{i});
    fprintf(' done \n')
end

for i=1:numel(htmlnames)
    fprintf('Conversion: mlx -> html (file %i of %i)...',i,numel(sourcefiles))
    matlab.internal.liveeditor.openAndConvert(sourcefiles{i},htmlfiles{i});
    fprintf(' done \n')
end

for i=1:numel(htmlnames)
    fprintf('Conversion: mlx -> m (file %i of %i)...',i,numel(sourcefiles))
    matlab.internal.liveeditor.openAndConvert(sourcefiles{i},mfiles{i});
    fprintf(' done \n')
end



