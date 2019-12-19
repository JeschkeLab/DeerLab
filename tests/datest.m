% datest    Test suite engine for DeerAnalysis
%
%   Usage:
%     datest            Run all tests
%     datest testname   Run all tests whose name starts with testname
%     
%   Options:
%       -d --display     Display graphical results of the tests
%       -r --regenerate  Recalculate and store regression data
%       -t --time        Report CPU timings of the tests
%       -p --perform     Report the maximal error found in the tests
%       -c --coverage    Run a coverage analysis over the tests
%       -b --badges      Generate JSON endpoint for badges
%       
%   Either the command syntax as above or the function syntax, e.g.
%   datest('asdf','t'), can be used. Any number of options are allowed
%
%   Run all tests including timings:   datest -t
%
%   All test files must have an underscore _ in their filename.

function out = datest(TestName,varargin)

% Check whether DeerAnalysis is on the Matlab path
DeerAnalysisPath = fileparts(which('datest'));
if isempty(DeerAnalysisPath)
    error('DeerAnalysis is not on the Matlab path!');
end

fid = 1; % output to command window

%Check for missing input arguments
if nargin<1
    TestName = '';
end
if nargin<2
    params = '';
end

%Accept more than one option as input
params = cell2mat(varargin);

if isempty(params)
    params = '';
end

%if not input test name is given but only options
if contains(['-d','-p','-r','-c','-t','-b'],TestName)
    params = [params TestName];
    TestName = '';
end
%Prepare options structure to be passed to the test functions
Opt.Display = contains(params,'-d') | contains(params,'--display');
Opt.Regenerate = contains(params,'-r') | contains(params,'--regen');
Opt.Verbosity = Opt.Display;

%Get internal options
displayErrors = contains(params,'-p') | contains(params,'--perform');
displayTimings = contains(params,'-t') | contains(params,'--time');
runCodeCoverage = contains(params,'-c') | contains(params,'--coverage');
makeBadges = contains(params,'-b') | contains(params,'--badge');

if any(TestName=='_')
    FileMask = [TestName '*.m'];
elseif strcmp(TestName,'all')
    FileMask = '*_*.m';
else
    FileMask = [TestName '*_*.m'];
end

%Look for test in the \tests directory of DeerAnalysis
FileList = dir(fullfile(DeerAnalysisPath,FileMask));

if numel(FileList)==0
    error('No test functions matching the pattern %s',FileMask);
end

TestFileNames = sort({FileList.name});

fprintf(fid,'=======================================================================\n');
fprintf(fid,'DeerAnalysis Unity Test Suite                 %s\n(Matlab %s)\n',datestr(now),version);
fprintf(fid,'DeerAnalysis location: %s\n',DeerAnalysisPath);
fprintf(fid,'=======================================================================\n');

% Codes for test outcomes:
%    0   test passed
%   +1   test failed
%   +2   test crashed
%   +3   not tested

OutcomeStrings = {'pass','failed','crashed','not tested'};

%get path to DeerAnalysis functions folder
path = fileparts(which('datest'));
path = path(1:end-length('\path'));

%list all the API functions, including models and private
Files1 = dir(fullfile(path,'functions','*.m'));
Files2 = dir(fullfile(path,'functions','models','*.m'));
Files = [Files1; Files2];
ExecutedLines = repmat({[]},length(Files),1);

%==========================================================================
%Loop over all tests to be run
%==========================================================================
for iTest = 1:numel(TestFileNames)
    
    if Opt.Display
        clf; drawnow;
    end
    
    thisTest = TestFileNames{iTest}(1:end-2);
    
    
    % Load, or regenerate, comparison data
    olddata = [];
    TestDataFile = ['data/' thisTest '.mat'];
    if exist(TestDataFile,'file')
        if Opt.Regenerate
            delete(TestDataFile);
            olddata = [];
        else
            try
                olddata = load(TestDataFile,'data');
                olddata = olddata.data;
            catch
                error('Could not load data for test ''%s''.',thisTest);
            end
        end
    end
    
    %Clear data in the cache before testing
    Pos = strfind(thisTest,'_');
    functionName = thisTest(1:Pos(1)-1);
    clear(functionName)
    
    % Clear and start profiler
    if runCodeCoverage
        profile clear
        profile on
    end
    
    tic
    try
        warning('off')
        if displayErrors
            [err,data,maxerr(iTest)] = feval(thisTest,Opt,olddata);
        else
            [err,data] = feval(thisTest,Opt,olddata);
        end
        warning('on')
        % if test returns empty err, then treat it as not tested
        if isempty(err)
            err = 3; % not tested
        else
            err = any(err~=0);
        end
        errorInfo = [];
        errorStr = '';
    catch exception
        data = [];
        err = 2;
        errorInfo = exception;
        errorStr = getReport(errorInfo);
        errorStr = ['    ' regexprep(errorStr,'\n','\n    ') char(10)];
    end
    time_used(iTest) = toc;
    
    %Retrieve profiler summary and turn profiler off
    if runCodeCoverage
        p = profile('info');
        profile off
        
        %Make list of all profiled function calls
        ExecutedFcns = [{p.FunctionTable(:).CompleteName}];
        %Analyze code coverage of each API function
        for n = 1:length(Files)
            FcnName = Files(n).name;
            pos = find(contains(ExecutedFcns,FcnName));
            if ~isempty(pos)
                %initialize containers
                for i=1:length(pos)
                    %get executed lines in profiler
                    tmp = p.FunctionTable(pos(i)).ExecutedLines;
                    container = ExecutedLines{n};
                    container(end+1:end+length(tmp(:, 1))) = tmp(:, 1);
                    ExecutedLines{n} = container;
                end
            end
        end
    end
    isRegressionTest = ~isempty(data);
    saveTestData = isRegressionTest && isempty(olddata);
    if saveTestData
        save(TestDataFile,'data');
    end
    
    testResults(iTest).err = double(err);
    testResults(iTest).err = double(err);
    testResults(iTest).name = thisTest;
    testResults(iTest).errorData = errorInfo;
    
    outcomeStr = OutcomeStrings{testResults(iTest).err+1};
    
    if ~isempty(data)
        typeStr = 'regression';
    else
        typeStr = 'direct';
    end
    
    if displayTimings
        timeStr = sprintf('%0.3f seconds',time_used(iTest));
    else
        timeStr = [];
    end
    
    try
        if displayErrors
            maxerrStr = sprintf('%0.2e ',maxerr(iTest));
        else
            maxerrStr = [];
        end
    catch
        maxerrStr = [];
    end
    str = sprintf('%-36s  %-12s%-8s%s%s\n%s',...
        testResults(iTest).name,typeStr,outcomeStr,maxerrStr,timeStr,errorStr);
    str(str=='\') = '/';
    
    testResults(iTest).msg = str;
    
    fprintf(fid,str);
    if Opt.Display
        if iTest<numel(TestFileNames), pause; end
    end
end

%Display results of code coverage analysis
if runCodeCoverage
    
    fprintf(fid,'-----------------------------------------------------------------------\n');
    fprintf(fid,'Code Coverage Analysis \n');
    fprintf(fid,'-----------------------------------------------------------------------\n');
    
    TotalCovered = 0;
    TotalRunnable = 0;
    %Analyze code coverage of each API function
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
        %account for end statement after return command
        if Runnable - Covered <= MissedEnds
            Covered = Runnable;
            Missed = [];
        end
        Coverage = 100*Covered/Runnable;
        %Print to console
        if (~isempty(TestName) && Coverage~=0) || isempty(TestName)
            if params =='c'
                fprintf('%-20s%-18s%3.2f%% %18s  %s \n',FcnName,' ',Coverage,'Lines missing:',mat2str(Missed))
            else
                fprintf('%-20s%-18s%3.2f%%\n',FcnName,' ',Coverage)
            end
        end
    end
    TotalCoverage = TotalCovered/TotalRunnable*100;
end
allErrors = [testResults.err];

% Display timings of slowest tests
if displayTimings
    fprintf(fid,'-----------------------------------------------------------------------\n');
    fprintf(fid,'Total test time:                        %7.3f seconds\n',sum(time_used));
    fprintf(fid,'Slowest tests:\n');
    [time,iTest] = sort(time_used,'descend');
    for q = 1:min(10,numel(time))
        fprintf(fid,'%-36s    %7.3f seconds\n',testResults(iTest(q)).name,time(q));
    end
end

% Display all tests that failed or crashed
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

msg = sprintf('%d passes, %d failures, %d crashes\n',Npasses,Nfails,Ncrashed);
fprintf(fid,msg);
if runCodeCoverage
    if isempty(TestName)
        fprintf('Total code coverage: %3.2f%%\n',TotalCoverage)
    end
end
fprintf(fid,'-----------------------------------------------------------------------\n');

% Return output if desired
if nargout==1
    out.Results = testResults;
    out.outcomes = allErrors;
end


if makeBadges
    
    %--------------------
    %Test Suite - Badge
    %--------------------
    
    %Determine grade of badge
    if Ncrashed~=0
        color = 'red';
    elseif Nfails~=0
        color = 'orange';
    else
        color = 'brightgreen';
    end
    %Write JSON endpoint file for shields.io
    fileID = fopen('testsuite_badge.json','w');
    fprintf(fileID,'{ \n');
    fprintf(fileID,'"schemaVersion": 1, \n');
    fprintf(fileID,'"label": "Tests", \n');
    fprintf(fileID,'"message": "%i pass, %i fail, %i crash", \n',sum(allErrors==0),sum(allErrors==1),sum(allErrors==2));
    fprintf(fileID,'"color": "%s"\n',color);
    fprintf(fileID,'} \n');
    fclose(fileID);
    
    if runCodeCoverage
        %--------------------
        %Code Coverage - Badge
        %--------------------
        %Determine grade of badge
        if TotalCoverage<50
            color = 'red';
        elseif TotalCoverage<75
            color = 'orange';
        else
            color = 'brightgreen';
        end
        %Write JSON endpoint file for shields.io
        fileID = fopen('coverage_badge.json','w');
        fprintf(fileID,'{ \n');
        fprintf(fileID,'"schemaVersion": 1, \n');
        fprintf(fileID,'"label": "Coverage", \n');
        fprintf(fileID,'"message": "%3.0f%%", \n',TotalCoverage);
        fprintf(fileID,'"color": "%s"\n',color);
        fprintf(fileID,'} \n');
        fclose(fileID);
    else
        warning('Coverage analysis was not requested. Coverage badge will no tbe generated.')
    end
end

return
