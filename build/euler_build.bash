#!/bin/bash

#Go to Euler root
cd ~    

#Prepare cloning of DeerAnalysis repo
echo Cloning DeerAnalysis repo...

#If there are any leftovers from a previous build, then remove them
if [ -d "~/dabuild" ] 
then
    echo "Removing old build files..."
	# rm -rf dabuild
fi
#Construct a fresh build directory
# mkdir -p dabuild
cd dabuild

#Clone repo
module load gcc/4.8.2 git/2.11.0 
# git clone https://github.com/luisfabib/DeerAnalysis2 
cd DeerAnalysis2
rm -rf .git
rm -rf .github
rm -rf .gitignore

echo "Launching MATLAB..."

module load new matlab/9.3

cd tests

matlab -singleCompThread -r CI_testing

exit