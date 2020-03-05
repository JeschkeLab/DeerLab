%  dltest    Test suite engine for DeerLab
% 
%    Usage:
%      dltest            Run all tests
%      dltest testname   Run all tests whose name starts with testname
%      
%    Options:
%        -d --display     Display graphical results of the tests
%        -r --regenerate  Recalculate and store regression data
%        -t --time        Report CPU timings of the tests
%        -p --perform     Report the maximal error found in the tests
%        -c --coverage    Run a coverage analysis over the tests
%        -b --badges      Generate JSON endpoint for badges
%        -u --tutorials   Run tutorial scripts
% 
%    Either the command syntax as above or the function syntax, e.g.
%    dltest('asdf','t'), can be used. Any number of options are allowed
% 
%    Run all tests including timings:   dltest -t
% 
%    All test files must have an underscore _ in their filename.

function out = dltest(TestName,varargin)

%  Check whether DeerLab is on the Matlab path
DeerLabPath = fileparts(which('dltest'));
if isempty(DeerLabPath)
    error('DeerLab is not on the Matlab path!');
end

fid = 1; %  output to command window

% Check for missing input arguments
if nargin<1
    TestName = '';
end
if nargin<2
    params = '';
end

% Accept more than one option as input
params = cell2mat(varargin);

if isempty(params)
    params = '';
end

% if not input test name is given but only options
if contains(TestName,'-')
    params = [params TestName];
    TestName = '';
end
% Prepare options structure to be passed to the test functions
Opt.Display = contains(params,'-d') | contains(params,'--display');
Opt.Regenerate = contains(params,'-r') | contains(params,'--regen');
Opt.Verbosity = Opt.Display;

% Get internal options
displayErrors = contains(params,'-p') | contains(params,'--perform');
displayTimings = contains(params,'-t') | contains(params,'--time');
runTutorials = contains(params,'-u') | contains(params,'--tutorials');
runCodeCoverage = contains(params,'-c') | contains(params,'--coverage');
makeBadges = contains(params,'-b') | contains(params,'--badge');

if any(TestName=='_')
    FileMask = [TestName '*.m'];
elseif strcmp(TestName,'all')
    FileMask = '*_*.m';
else
    FileMask = [TestName '*_*.m'];
end

% Look for test in the \tests directory of DeerLab
FileList = dir(fullfile(DeerLabPath,FileMask));

if numel(FileList)==0
    error('No test functions matching the pattern % s',FileMask);
end

TestFileNames = sort({FileList.name});

fprintf(fid,'=======================================================================\n');
fprintf(fid,'DeerLab Unit Test Suite                 % s\n(Matlab % s)\n',datestr(now),version);
fprintf(fid,'DeerLab location: % s\n',DeerLabPath);
fprintf(fid,'=======================================================================\n');

%  Codes for test outcomes:
%     0   test passed
%    +1   test failed
%    +2   test crashed
%    +3   not tested

OutcomeStrings = {'pass','failed','crashed','not tested'};

% get path to DeerLab functions folder
path = fileparts(which('dltest'));
path = path(1:end-length('\path'));

% list all the API functions, including models and private
Files1 = dir(fullfile(path,'functions','*.m'));
Files2 = dir(fullfile(path,'functions','models','*.m'));
Files = [Files1; Files2];
ExecutedLines = repmat({[]},length(Files),1);

% ==========================================================================
% Loop over all tests to be run
% ==========================================================================
for iTest = 1:numel(TestFileNames)
    
    if Opt.Display
        clf; drawnow;
    end
    
    thisTest = TestFileNames{iTest}(1:end-2);
    
    % Clear data in the cache before testing
    Pos = strfind(thisTest,'_');
    functionName = thisTest(1:Pos(1)-1);
    clear(functionName)
    
    %  Clear and start profiler
    if runCodeCoverage
        profile clear
        profile on
    end
    
    % Run test, catch any errors
    testFcn = str2func(thisTest);
    nArgsOut = nargout(testFcn);
    nArgsIn = nargin(testFcn);
    warning('off')
    tic
    try
         
        maxerr(iTest) = 0;
        if nArgsOut==1
            if nArgsIn==0
                pass = testFcn();
            else
                pass = testFcn(Opt);
            end
        else
            if nArgsIn<1
                error('1 input is needed.');
            end
            if nArgsOut==2
                [pass,maxerr(iTest)] = testFcn(Opt);
            else
                [pass] = testFcn(Opt);
            end
        end
        % If test returns empty err, then treat it as not tested
        if isempty(pass)
            pass = 3; %  not tested
        else
            pass = any(pass~=1);
        end
        errorInfo = [];
        errorStr = '';
    catch exception
         
        pass = 2;
        errorInfo = exception;
        errorStr = getReport(errorInfo);
        errorStr = ['    ' regexprep(errorStr,'\n','\n    ') newline];
    end
    time_used(iTest) = toc;
    warning('on')
    
    % Retrieve profiler summary and turn profiler off
    if runCodeCoverage
        p = profile('info');
        profile off
        
        % Make list of all profiled function calls
        ExecutedFcns = [{p.FunctionTable(:).CompleteName}];
        % Analyze code coverage of each API function
        for n = 1:length(Files)
            FcnName = Files(n).name;
            pos = find(contains(ExecutedFcns,FcnName));
            if ~isempty(pos)
                % initialize containers
                for i=1:length(pos)
                    % get executed lines in profiler
                    tmp = p.FunctionTable(pos(i)).ExecutedLines;
                    container = ExecutedLines{n};
                    container(end+1:end+length(tmp(:, 1))) = tmp(:, 1);
                    ExecutedLines{n} = container;
                end
            end
        end
    end
    
    testResults(iTest).err = double(pass);
    testResults(iTest).err = double(pass);
    testResults(iTest).name = thisTest;
    testResults(iTest).errorData = errorInfo;
    
    outcomeStr = OutcomeStrings{testResults(iTest).err+1};
    
    if displayTimings
        timeStr = sprintf('% 0.3f sec',time_used(iTest));
    else
        timeStr = [];
    end
    
    try
        if displayErrors
            if ~isnan(maxerr(iTest))
                maxerrStr = sprintf('% 0.2e ',maxerr(iTest));
            else
                maxerrStr = ' --------';
            end
        else
            maxerrStr = [];
        end
    catch
        maxerrStr = [];
    end
   
    str = sprintf('% -38s   %-8s% -8s% s \n',...
        testResults(iTest).name,outcomeStr,maxerrStr,timeStr);
    str(str=='\') = '/';
    
    
    str2 = sprintf('% -38s   %-8s% -8s% s \n% s',...
        testResults(iTest).name,outcomeStr,maxerrStr,timeStr,errorStr);
    str2(str2=='\') = '/';
    testResults(iTest).msg = str2;
    
    switch outcomeStr
        case 'pass'
            fid = 1;
        case 'failed'
            fid = 1;
            str = ['[\b' str ']\b'];
        case 'crashed'
            fid = 2;
    end
    
    fprintf(fid,str);
    
    if Opt.Display
        if iTest<numel(TestFileNames), pause; end
    end
end

fid = 1;

allErrors = [];

% Test that tutorials run without errors
if runTutorials
    
    fprintf(fid,'-----------------------------------------------------------------------\n');
    fprintf(fid,'Tutorials tests \n');
    fprintf(fid,'-----------------------------------------------------------------------\n');
    
    % Get contents of tutorials directory
    files = dir(fullfile(erase(DeerLabPath,'tests'),'/tutorials/**/*.mlx'));
    sourcenames = {files.name};
    mnames = cellfun(@(S)[S(1:end-3) 'm'],sourcenames,'UniformOutput',false);
    
    % Generate full paths
    paths = {files.folder};
    sourcefiles = cellfun(@(x,y)fullfile(x,y),paths,sourcenames,'UniformOutput',false);
    path = fullfile(DeerLabPath,'_tutorials');
    if ~exist(path,'dir')
        mkdir(path)
    end
    addpath(genpath(fullfile(erase(DeerLabPath,'tests'),'/tutorials')));
    mfiles = cellfun(@(y)fullfile(path,y),mnames,'UniformOutput',false);
    
    for i=1:numel(mfiles)
        % If m-files are not present then compile the mlx files in dummy directory
        matlab.internal.liveeditor.openAndConvert(sourcefiles{i},mfiles{i});
        % Read the m-file
        fid2 = fopen(mfiles{i},'r');   
        f = fread(fid2,'*char')';
        fclose(fid2);
        % If there are any clear, clc or clf remove them
        f = strrep(f,'clear','%');
        f = strrep(f,'clc','%');
        f = strrep(f,'clf','%');
        % Save the modified m-file
        fid2 = fopen(mfiles{i},'w+');   
        fprintf(fid2,'%s',f);
        fclose(fid2);
    end
    
    % Change into dummy directory
    currentdir = pwd;
    cd(fileparts(mfiles{1}));
    % Disable figure/plot generation during tests
    set(0,'DefaultFigureVisible','off')

    % Run all tutorial scripts
    for i=1:numel(mfiles)
        % Skip the template file
        if contains(mfiles{i},'template.m')
            continue
        end
        % Run turorial, catch errors
        try
            tic
            % Run script without print outputs
            evalc(erase(mnames{i},'.m'));
            outcomeStr = 'pass';
            errorStr = '';
            allErrors(end+1) = 0;
            fid = 1;
        catch exception
            outcomeStr = 'crash';
            allErrors(end+1) = 2;
            errorStr = getReport(exception);
            errorStr = ['    ' regexprep(errorStr,'\n','\n    '),newline];
            fid = 2;
        end
        tutorialTime = toc;
        
        if displayTimings
            timeStr = sprintf('% 0.3f seconds',tutorialTime);
        else
            timeStr = [];
        end
        
        % Print result
        str = sprintf('% -36s  % -8s %s \n %s',mnames{i},outcomeStr,timeStr,errorStr);
        str(str=='\') = '/';
        fprintf(fid,str)
    end
    % Restore figure/plot generation
    set(0,'DefaultFigureVisible','on')
    %Get back into original directory
    cd(currentdir);
    
end

% Display results of code coverage analysis
if runCodeCoverage
    
    fprintf(fid,'-----------------------------------------------------------------------\n');
    fprintf(fid,'Code Coverage Analysis \n');
    fprintf(fid,'-----------------------------------------------------------------------\n');
    
    TotalCovered = 0;
    TotalRunnable = 0;
    % Analyze code coverage of each API function
    for n = 1:length(Files)
        FcnName = Files(n).name;
        Path = Files(n).folder;
        RunnableLines = callstats('file_lines',fullfile(Path,FcnName));
        TotalRunnable = TotalRunnable + length(unique(RunnableLines));
        Executed = unique(ExecutedLines{n});
        Covered = length(Executed);
        TotalCovered = TotalCovered + Covered;
        Runnable = length(unique(RunnableLines));
        Code = fileread(FcnName);
        if params =='l'
            Missed = RunnableLines;
            for k=1:length(Executed)
                Missed(RunnableLines==Executed(k)) = NaN;
            end
            Missed(isnan(Missed)) = [];
        end
        MissedEnds = length(strfind(Code,'error')) + length(strfind(Code,'return')) ...
            + length(strfind(Code,'break')) + length(strfind(Code,'fprintf')) + length(strfind(Code,'plot'));
        % account for end statement after return command
        if Runnable - Covered <= MissedEnds
            Covered = Runnable;
            Missed = [];
        end
        Coverage = 100*Covered/Runnable;
        % Print to console
        if (~isempty(TestName) && Coverage~=0) || isempty(TestName)
            if params =='c'
                fprintf('%-20s%-18s%6.2f%18s    %s\n',FcnName,' ',Coverage,'Lines missing:',mat2str(Missed))
            else
                fprintf('%-20s%-18s%6.2f\n',FcnName,' ',Coverage);
            end
        end
    end
    TotalCoverage = TotalCovered/TotalRunnable*100;
end
allErrors = [allErrors [testResults.err]];

%  Display timings of slowest tests
if displayTimings
    fprintf(fid,'-----------------------------------------------------------------------\n');
    fprintf(fid,'Total test time:                        % 7.3f seconds\n',sum(time_used));
    fprintf(fid,'Slowest tests:\n');
    [time,iTest] = sort(time_used,'descend');
    for q = 1:min(10,numel(time))
        fprintf(fid,'% -36s    % 7.3f seconds\n',testResults(iTest(q)).name,time(q));
    end
end

%  Display all tests that failed or crashed
if any(allErrors==1) || any(allErrors==2)
    fprintf(fid,'-----------------------------------------------------------------------\n');
    for iTest = find(allErrors)
        fprintf(fid,testResults(iTest).msg);
    end
end

fprintf(fid,'-----------------------------------------------------------------------\n');

Ncrashed = sum(allErrors==2);
Npasses = sum(allErrors==0);
Nfails = sum(allErrors==1);

msg = sprintf('% d passes, % d failures, % d crashes\n',Npasses,Nfails,Ncrashed);
fprintf(fid,msg);
if runCodeCoverage
    if isempty(TestName)
        fprintf('Total code coverage: % 3.2f% % \n',TotalCoverage)
    end
end
fprintf(fid,'-----------------------------------------------------------------------\n');

%  Return output if desired
if nargout==1
    out.Results = testResults;
    out.Outcomes = allErrors;
    out.Errors = sum(allErrors);
end


if makeBadges
    
    % --------------------
    % Test Suite - Badge
    % --------------------
    
    % Determine grade of badge
    if Ncrashed~=0
        color = 'red';
    elseif Nfails~=0
        color = 'orange';
    else
        color = 'brightgreen';
    end
    % Write JSON endpoint file for shields.io
    fileID = fopen('testsuite_badge.json','w');
    fprintf(fileID,'{ \n');
    fprintf(fileID,'"schemaVersion": 1, \n');
    fprintf(fileID,'"label": "Tests", \n');
    fprintf(fileID,'"message": "% i pass, % i fail, % i crash", \n',sum(allErrors==0),sum(allErrors==1),sum(allErrors==2));
    fprintf(fileID,'"color": "% s"\n',color);
    fprintf(fileID,'} \n');
    fclose(fileID);
    fprintf('Test suite badge created: % s\n',fullfile(pwd,'testsuite_badge.json'))
    
    if runCodeCoverage
        % --------------------
        % Code Coverage - Badge
        % --------------------
        % Determine grade of badge
        if TotalCoverage<50
            color = 'red';
        elseif TotalCoverage<75
            color = 'orange';
        else
            color = 'brightgreen';
        end
        % Write JSON endpoint file for shields.io
        fileID = fopen('coverage_badge.json','w');
        fprintf(fileID,'{ \n');
        fprintf(fileID,'"schemaVersion": 1, \n');
        fprintf(fileID,'"label": "Coverage", \n');
        fprintf(fileID,'"message": "%3.0f %% ", \n',TotalCoverage);
        fprintf(fileID,'"color": "%s"\n',color);
        fprintf(fileID,'} \n');
        fclose(fileID);
        fprintf('Coverage badge created: % s\n',fullfile(pwd,'coverage_badge.json'))
    else
        warning('Coverage analysis was not requested. Coverage badge will no tbe generated.')
    end
    
    fprintf(fid,'-----------------------------------------------------------------------\n');
    
end

return
