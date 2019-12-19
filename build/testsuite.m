
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
!python uploadS3.py --file "coverage_badge.json" --bucket deershields
!python uploadS3.py --file "testsuite_badge.json" --bucket deershields

%Stop diary function
diary off