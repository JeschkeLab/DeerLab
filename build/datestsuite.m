
%Stop diary function in case it is running
diary off

%Add paths to dependencies
addpath('../functions')
addpath('../functions/models')
addpath('../tests')
addpath('../tests/comparison')

%Log all console outputs to artifact log
logname = [datestr(now,'yyyymmdd_HHMMSS'),'_testsuite_log'];
diary(logname)
diary on

%Run full test suite 
datest --time --perform --coverage --badges

%Load the badges JSON endpoints to AWS
if ispc
    system('python uploadS3.py --keyfile %USERPROFILE%\.ssh\aws_access_keys.txt --file "coverage_badge.json" --bucket deershields')
    system('python uploadS3.py --keyfile %USERPROFILE%\.ssh\aws_access_keys.txt --file "testsuite_badge.json" --bucket deershields')
elseif isunix
    %system('python3 uploadS3.py --keyfile ~/.ssh/aws_access_keys.txt --file "coverage_badge.json" --bucket deershields')
    %system('python3 uploadS3.py --keyfile ~/.ssh/aws_access_keys.txt --file "testsuite_badge.json" --bucket deershields')
elseif ismac
    error('macOS is not supported')
end

%Remove artifacts once uploaded
% delete coverage_badge.json
% delete testsuite_badge.json

%Stop diary function
diary off