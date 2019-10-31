import glob
import os
import datetime
import shutil
import boto3
import boto
import time
import pytz
from time import mktime

def getListOfFiles(dirName):
    #Create a list of file and sub directories 
    #Names in the given directory 
    listOfFile = os.listdir(dirName)
    allFiles = list()
    #Iterate over all the entries
    for entry in listOfFile:
        #Create full path
        fullPath = '/'.join([dirName, entry])
        #If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
                
    return allFiles

#Load and read the secure access keys to connect to the AWS S3 client
if os.path.isfile(".\aws_access_keys.txt"):
	print("AWS access keys not found. You may not have rights to request this action.")
	exit()

with open('aws_access_keys.txt', 'r') as file:
    AccessKeys = [line.rstrip('\n') for line in file]

# Create an S3 client
print("Establishing AWS S3 client connection")
s3_client = boto3.client('s3',
		 aws_access_key_id = AccessKeys[0],
		 aws_secret_access_key = AccessKeys[1])
		 
#Get the keys to the DeerAnalysis bucket
print("Retrieving S3 bucket information")
s3 = boto3.resource('s3',
		 aws_access_key_id = AccessKeys[0],
		 aws_secret_access_key = AccessKeys[1],
		 region_name='eu-west-1')
bucket = s3.Bucket('deeranalysis.org')

base_dir = "../docs"
#Get full list of local files
localFiles = getListOfFiles(base_dir)

for key in bucket.objects.all():
		
		#Get full local path of current file in bucket 
		filepath = '/'.join([base_dir, key.key])
		
		#Remove bucket file from list of local files
		if filepath in localFiles: localFiles.remove(filepath)
		
		#If file has been removed from the local source, remove it from the bucket
		if not os.path.isfile(filepath):
			print("Removing ",filepath," from bucket, not found in local source...")
			s3_client.delete_object(Bucket='deeranalysis.org',Key=key.key)
			continue
		
		#Get last modified date of local files
		modifyDate = datetime.datetime.fromtimestamp(os.path.getmtime(filepath))
		
		#Localize timestamp of local files to West-Europe timezone
		euwest = pytz.timezone('Europe/Amsterdam')
		modifyDate = euwest.localize(modifyDate)
		
		filePathList = filepath.split("/") 
		filename = filePathList[-1] #The last element is a the filename
		
		#Update the file if the local source file is newer than the version in the S3 bucket
		if modifyDate > key.last_modified:
			print("Updating", filepath, " in web bucket... ")
			s3_client.upload_file(filepath,'deeranalysis.org',filename)

#Add the remaining local files which are still not on the we bucket
for files in localFiles:
	filename = files.split("/")
	filename = filename[-1] 
	print("Adding", filepath, " to web bucket... ")
	s3_client.upload_file(files,'deeranalysis.org',filename)	

print("Finished: AWS S3 DeerAnalysis.org bucket is up to date.")