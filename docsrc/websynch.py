import glob
import os
import datetime
import shutil
import boto3
import time
import pytz
from time import mktime
from argparse import ArgumentParser
from pathlib import Path
 
#Parse options if access keys are being passed as inputs
parser = ArgumentParser()
parser.add_argument("-k", "--key",dest="accesskey",default=[],help="AWS Access Key")
parser.add_argument("-s", "--secretkey",dest="secretkey", default=[],help="AWS Secret Access Key")
args = parser.parse_args()

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

def metadataType(file):
		if file.endswith('.html'):
				metastr =  "text/html"
		elif file.endswith('.css'):
				metastr =  "text/css"
		elif file.endswith('.js'):
				metastr =  "application/javascript"
		elif file.endswith('.svg'):
				metastr =   "image/svg+xml"
		elif file.endswith('.png'):
				metastr =  "image/png"
		elif file.endswith('.gif'):
				metastr =  "image/gif"
		elif file.endswith('.jpg'):
				metastr =  "image/jpg"
		else: 
				metastr = ""
		
		return metastr


#If access codes have not been passed as inputs then check if file is there
if not args.accesskey or not args.secretkey:

	#Load and read the secure access keys to connect to the AWS S3 client
	if os.path.isfile(".\aws_access_keys.txt"):
		print("AWS access keys not found. You may not have rights to request this action.")
		exit()

	with open('aws_access_keys.txt', 'r') as file:
		AccessKeys = [line.rstrip('\n') for line in file]
else:
	 AccessKeys = [args.accesskey, args.secretkey]

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
		file = '/'.join([base_dir, key.key])
		
		#Remove bucket file from list of local files
		if file in localFiles: localFiles.remove(file)
		
		#If file has been removed from the local source, remove it from the bucket
		if not os.path.isfile(file):
			print("Removing ",file," from bucket, not found in local source...")
			s3_client.delete_object(Bucket='deeranalysis.org',Key = key.key)
			continue
		
		#Get last modified date of local files
		modifyDate = datetime.datetime.fromtimestamp(os.path.getmtime(file))
		
		#Localize timestamp of local files to West-Europe timezone
		euwest = pytz.timezone('Europe/Amsterdam')
		modifyDate = euwest.localize(modifyDate)
		
		#Update the file if the local source file is newer than the version in the S3 bucket
		if modifyDate > key.last_modified:
			print("Updating", file, "in web bucket... ")
			s3.meta.client.upload_file(file, 'deeranalysis.org', key.key, ExtraArgs={'ContentType': metadataType(file)} )

#Add the remaining local files which are still not on the we bucket
for files in localFiles:
	key = files.replace("../docs/","")
	print("Adding", key, "to web bucket... ")
	s3.meta.client.upload_file(file, 'deeranalysis.org', key.key, ExtraArgs={'ContentType': metadataType(file)} )

print("Finished: AWS S3 DeerAnalysis.org bucket is up to date.")