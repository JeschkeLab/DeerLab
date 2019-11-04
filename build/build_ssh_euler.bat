
@ECHO OFF


echo.Uploading DeerAnalysis source to Euler...
scp -q -r ../../DeerAnalysis/docsrc luisf@euler.ethz.ch:~/DeerAnalysis/docsrc
scp -q -r ../../DeerAnalysis/functions luisf@euler.ethz.ch:~/DeerAnalysis/functions
scp -q -r ../../DeerAnalysis/tests luisf@euler.ethz.ch:~/DeerAnalysis/tests
scp -q -r ../../DeerAnalysis/build luisf@euler.ethz.ch:~/DeerAnalysis/build
echo.All files required for build have been successfully uploaded.

ssh luisf@euler.ethz.ch `~/DeerAnalysis/build/euler_build.bash`

