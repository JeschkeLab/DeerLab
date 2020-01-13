import glob
import os
import datetime
import shutil
import boto3
import time
import pytz
import pathlib
from time import mktime
from argparse import ArgumentParser
from pathlib import Path
 
#Parse options if access keys are being passed as inputs
parser = ArgumentParser()
parser.add_argument("-i", "--keyfile",dest="keyfile",default=[],help="AWS Access Key")
parser.add_argument("-k", "--key",dest="accesskey",default=[],help="AWS Access Key")
parser.add_argument("-s", "--secretkey",dest="secretkey", default=[],help="AWS Secret Access Key")
parser.add_argument("-f", "--file",dest="filename", default=[],help="File to upload")
parser.add_argument("-b", "--bucket",dest="bucketname", default=[],help="Name of S3 bucket")
args = parser.parse_args()

def metadataType(file):
		if file.endswith('.html'):
				metastr =  "text/html"
		elif file.endswith('.css'):
				metastr =  "text/css"
		elif file.endswith('.js'):
				metastr =  "application/javascript"
		elif file.endswith('.json'):
				metastr =  "application/json"
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
	
	if not args.keyfile:
		keypath = 'aws_access_keys.txt'
	else:
		keypath = args.keyfile
	
	#Set absolute path to key file
	#cwd = os.getcwd()
	#keypath =  os.path.join(cwd,keypath)

	#Load and read the secure access keys to connect to the AWS S3 client
	if not os.path.isfile(keypath):
		print("AWS access keys not found. You may not have rights to request this action.")
		exit()
	with open(keypath, 'r') as file:
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
bucket = s3.Bucket(args.bucketname)

#Get filename to upload
localFile = args.filename

print("Uploading ", localFile, "to S3 bucket... ")
s3.meta.client.upload_file(localFile, args.bucketname, localFile, ExtraArgs={'ContentType': "application/json"} )

print("Finished: Data transfer successful.")