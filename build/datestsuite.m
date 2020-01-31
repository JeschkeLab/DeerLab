%Add paths to dependencies
addpath('../functions')
addpath('../functions/models')
addpath('../tests')
addpath('../tests/comparison')

%Run full test suite 
Out = datest('--time','--perform','--coverage','--badges');

if Out.Errors>0
    logname = 'datestsuite.error';
    save(logname)
end

%Load the badges JSON endpoints to AWS
if ispc
    system('python uploadS3.py --keyfile %USERPROFILE%\.ssh\aws_access_keys.txt --file "coverage_badge.json" --bucket deershields')
    system('python uploadS3.py --keyfile %USERPROFILE%\.ssh\aws_access_keys.txt --file "testsuite_badge.json" --bucket deershields')

    %Delete files once uploaded
    delete('coverage_badge.json')
    delete('testsuite_badge.json')

elseif isunix
    %system('python3 uploadS3.py --keyfile ~/.ssh/aws_access_keys.txt --file "coverage_badge.json" --bucket deershields')
    %system('python3 uploadS3.py --keyfile ~/.ssh/aws_access_keys.txt --file "testsuite_badge.json" --bucket deershields')
elseif ismac
    error('macOS is not supported')
end