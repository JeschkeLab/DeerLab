function docs(pagename)

DApath = fileparts(which('docs'));
DApath = DApath(1:end-numel('\functions'));
HTMLfile = fullfile(DApath,'docs\api',[pagename '.html']);
HTMLfile2 = fullfile(DApath,'docs\models',[pagename '.html']);

%Check whether file exists
if exist(HTMLfile,'file')
    %Then open the corresponding HMLT file using the MATLAB browser
    web(HTMLfile)
elseif exist(HTMLfile2,'file')
    %Then open the corresponding HMLT file using the MATLAB browser
    web(HTMLfile2)
else
    %get path to DeerAnalysis functions folder
    path = fileparts(which('datest'));
    path = path(1:end-length('\path'));
    %list all the API functions, including models and private
    Files1 = dir(fullfile(path,'functions','*.m'));
    Files2 = dir(fullfile(path,'functions','models','*.m'));
    Files = [Files1; Files2];
    FcnName = {Files(:).name};
    
    %Look for close matches and recommend to the user
    pos = [];
    charmatch = 5;
    while isempty(pos) && charmatch>3
        nearbynames = strncmpi(strtrim(pagename),FcnName,charmatch);
        pos = find(nearbynames);
        charmatch = charmatch - 1;
    end
    if ~isempty(pos)
        fprintf('Function name not found. Did you mean one of these? \n')
        for i=1:length(pos)
            name = FcnName{pos(i)};
            name = name(1:end-2);
            fprintf('%20s \n',name)
        end
    else
        fprintf('Function name not found. No close matches could be found. \n')
    end
    
end