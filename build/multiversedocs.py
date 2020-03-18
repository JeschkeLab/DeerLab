import subprocess
import fileinput
import glob
import os
import shutil

def formatProcOut(input):
    output = str(input)
    output = output.replace("b'","")
    output = output.replace("'","")
    return output

ignoredtags = ['v.0.1-beta','0.1-beta','0.2.beta','0.5.beta']

#Get list of all tags in the DeerLab repo
tags = subprocess.check_output (["git", "tag"])
#Format output 
tags = formatProcOut(tags)
tags = tags.split("\\n")
tags = tags[:-2]
tags.insert(0,"develop")

for tag in tags: 
    if any(tag in ignore for ignore in ignoredtags):
        tags.remove(tag)

path = '../multidocs'
if os.path.exists(path):
    shutil.rmtree(path)
    
#Get tag of most modern version
mostrecent = subprocess.check_output (["git", "describe", "--tags", "--abbrev=0"])
#Format output
mostrecent = formatProcOut(mostrecent)
mostrecent = mostrecent.replace("\\n","")

startbranch = subprocess.check_output (["git", "branch", "--show-current"])
#Format output
startbranch = formatProcOut(startbranch)
startbranch = startbranch.replace("\\n","")

for tag in tags:

    if tag == 'develop':
        subprocess.run (["git", "checkout", "develop"])
        #Build source code in .\docsrc
        subprocess.run (["..\docsrc\make.bat", "clean"])
    else:
        subprocess.run (["git", "checkout", "master"])
        #Get commit SHA corresponding to current tag
        commit = subprocess.check_output (["git", "rev-list", "-n", "1", tag])
        commit = formatProcOut(commit)
        commit = commit.replace("\\n","")
        #Checkout that commit
        subprocess.run (["git", "checkout","-f", commit])
        #Build source code in .\docsrc
        subprocess.run (["..\docsrc\make.bat"])

    #The devleopment version is compiled first, then copy the development index.html to the rest
    if tag != 'develop':
         shutil.copyfile('../multidocs/develop/index.html', '../docs/index.html')
         counter = 0
         for line in fileinput.input('../docs/index.html', inplace=True):
            if not counter:
                if '<div class="select">' in line:
                    counter = 1000
                else:
                    print(line,end = '')
            else:
                if '<div class="select_arrow">\n' in line:
                    counter = 3
                counter -= 1
    #Get list from all HTML files in the compiled documentation
    path = '../docs/'
    htmlfiles = [f for f in glob.glob(path + "**/*.html", recursive = True)]


    #Loop over all HTML available at that version
    for file in htmlfiles:
    

        relpath = file.replace('../docs','.')
        
        print('Processing ' + file + '...                                     ', end='\r')
        #Go through all files
        for line in fileinput.FileInput(file,inplace=1):
                    
            #Find the place to put the HTML code for the version dropdown menu
            if '<div role="search">' in line:
                line = line.replace(line,'</div>\n\n\n' + line)
                line = line.replace(line,'    </div>\n' + line)
                line = line.replace(line,'    <div class="select_arrow">\n' + line)
                line = line.replace(line,'    </select>\n' + line)
                #Make the current version the default selection
                for version in tags:
                            
                     
                    href = os.path.relpath('multidocs/' + version + relpath, 'multidocs/' + tag + relpath)
                    href = href[3:]
                    if tag == version:
                        line = line.replace(line,'        <option selected value="' + href + '">' + version + '</option>\n' + line)
                    else:
                        line = line.replace(line,'        <option value="' + href + '">' + version + '</option>\n' + line)
                line = line.replace(line,'    <select onChange="window.location.href=this.value">\n' + line)
                line = line.replace(line,'<div class="select">\n' + line)
            
            if 'custom.css' in line:
                if tag != 'develop':
                    href = os.path.relpath('multidocs/' + tag + relpath, 'multidocs/' + 'develop' + relpath)
                    href = href[3:]
                    idx = href.find(tag)
                    href = href[:idx] + 'develop/_static/custom.css'
                    line = line.replace(line,'<link rel="stylesheet" href="' + href + '" type="text/css" />')
            
            if 'theme_override.css' in line and tag != 'develop':
                    href = os.path.relpath('multidocs/' + tag + relpath, 'multidocs/' + 'develop' + relpath)
                    href = href[3:]
                    idx = href.find(tag)
                    href = href[:idx] + 'develop/_static/theme_override.css'
                    line = line.replace(line,'<link rel="stylesheet" href="' + href + '" type="text/css" />')   
                    
            # Edits in the relesed versions index.html files
            if '<title>DeerLab development-version documentation' in line and tag != 'develop':
                    line = line.replace(line,'<title>DeerLab documentation</title>')
            if '<li>DeerLab development-version documentation' in line and tag != 'develop':
                    line = line.replace(line,'<li>DeerLab '+tag+' documentation</li>')    

            if '<h1>DeerLab development-version documentation' in line and tag != 'develop':
                    line = line.replace(line,'<h1>DeerLab '+tag+' documentation</h1>')  
                    
            #Print line into file
            print(line,end = '')
            

    #Copy current build into collective docs
    if tag == mostrecent:
        shutil.copytree('../docs', '../multidocs/' + tag)
    else:
        shutil.copytree('../docs', '../multidocs/' + tag)
    
    #Copy the homepage html file to the multidocs build
    if tag == 'develop':
        shutil.copyfile('../docsrc/source/homepage/index.html', '../multidocs/index.html')
        shutil.copytree('../docsrc/source/homepage/_static', '../multidocs/_static')
        
    print('File processing for version ' + tag + ' completed                                \n', end='\r')
    print('\n')



#Return to the branch the script started at
subprocess.run (["git", "checkout", startbranch])


