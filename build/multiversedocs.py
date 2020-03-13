import subprocess
import fileinput
import glob
import os

def formatProcOut(input):
	output = str(input)
	output = output.replace("b'","")
	output = output.replace("'","")
	return output


#Get list from all HTML files in the compiled documentation
path = '../docs/'
htmlfiles = [f for f in glob.glob(path + "**/*.html", recursive=True)]

#Get list of all tags in the DeerLab repo
tags = subprocess.check_output (["git", "tag"])
#Format the output 
tags = formatProcOut(tags)
tags = tags.split("\\n")
tags = tags[:-2]

for tag in tags:

	os.system('..\docsrc\make.bat clean')

	commit = subprocess.check_output (["git", "rev-list", "-n", "1", tag])
	
	#Get tag of current version
	currentag = subprocess.check_output (["git", "describe", "--tags", "--abbrev=0"])
	#Format output
	currentag = formatProcOut(currentag)
	currentag = currentag.replace("\\n'","")

	#Loop over all HTML available at that version
	for file in htmlfiles:
	
		print('Processing '+file+'...                                     ', end='\r')
		#Go through all files
		for line in fileinput.FileInput(file,inplace=1):
		
			#Find the place to put the code for the version dropdown menu
			if '<div role="search">' in line:
				line=line.replace(line,'</div>\n\n\n' + line)
				line=line.replace(line,'    </div>\n' + line)
				line=line.replace(line,'    <div class="select_arrow">\n' + line)
				line=line.replace(line,'    </select>\n' + line)
				#Make the current version the default selection
				for version in tags:
					if currentag == version:
						line=line.replace(line,'        <option selected>'+version+'</option>\n' + line)
					else:
						line=line.replace(line,'        <option>'+version+'</option>\n' + line)
				line=line.replace(line,'    <select>\n' + line)
				line=line.replace(line,'<div class="select">\n' + line)
			#Print lines into file
			print(line,end='')
			

	print('File processing for version '+tag+' completed                                \n', end='\r')
