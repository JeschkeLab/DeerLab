import subprocess
import fileinput
import glob
import os
import shutil
import sys

def formatProcOut(input):
    output = str(input)
    output = output.replace("b'","")
    output = output.replace("'","")
    return output

ignoredtags = ['v.0.1-beta','0.1-beta','0.2.beta','0.5.beta','0.6.beta','0.8.beta']

#Get list of all tags in the DeerLab repo
tags = subprocess.check_output (["git", "tag"])
#Format output 
tags = formatProcOut(tags)
tags = tags.split("\\n")
tags.insert(0,"develop")

newtags = list()
for tag in tags: 
    if not any(tag in ignore for ignore in ignoredtags):
        newtags.append(tag)
tags = newtags;

path = os.path.join('..','multidocs')
if os.path.exists(path):
    shutil.rmtree(path)

#Get tag of most modern version
mostrecent = subprocess.check_output (["git", "describe", "--tags", "--abbrev=0"])
#Format output
mostrecent = formatProcOut(mostrecent)
mostrecent = mostrecent.replace("\\n","")

startbranch = subprocess.check_output (["git", "branch", "--show-current"]) # This requires git 2.2 or newer
#Format output
startbranch = formatProcOut(startbranch)
startbranch = startbranch.replace("\\n","")

for tag in tags:

    cddir = os.path.join('..','docsrc')
    os.chdir (cddir)
    src = os.path.join('.', 'source')
    dest = os.path.join('..', 'docs')
    #if os.path.exists(dest):
    #    shutil.rmtree(dest)
    
    buildcmd = ['sphinx-build', '-Q','-E', '-b', 'html',src,dest]
    if tag == 'develop':
        subprocess.run (["git", "checkout", "develop"])
        #Build source code in .\docsrc
        subprocess.run (buildcmd)
    else:
        subprocess.run (["git", "checkout", "master"])
        #Get commit SHA corresponding to current tag
        commit = subprocess.check_output (["git", "rev-list", "-n", "1", tag])
        commit = formatProcOut(commit)
        commit = commit.replace("\\n","")
        #Checkout that commit
        subprocess.run (["git", "checkout","-f", commit])
        #Build source code in .\docsrc
        subprocess.run (buildcmd)
        
    cddir = os.path.join('..','build')
    os.chdir (cddir)

    #The devleopment version is compiled first, then copy the development index.html to the rest
    if tag != 'develop':
        src = os.path.join('..','multidocs','develop','index.html')
        dest = os.path.join('..','docs','index.html')
        shutil.copyfile(src, dest)
        counter = 0
        path = os.path.join('..','docs','index.html')
        for line in fileinput.input(path, inplace=True,mode='rb'):
            
            line = line.decode('utf-8')
            if not counter:
                if '<div class="select">' in line:
                    counter = 1000
                else:
                    sys.stdout.buffer.write(line.encode('utf-8'))
            else:
                if '<div class="select_arrow">\n' in line:
                    counter = 3
                counter -= 1
    
    #Get list from all HTML files in the compiled documentation
    docspath = os.path.join('..','docs')
    multidocspath = os.path.join('..','multidocs')
    filterpath = os.path.join('**','*.html')
    htmlfiles = [f for f in glob.glob(os.path.join(docspath,filterpath), recursive = True)]
    
    #Loop over all HTML available at that version
    for file in htmlfiles:
        
        relpath = file.replace(docspath,'.')

        print('Processing ' + file + '...                                     ', end='\r')
        
        #Go through all lines
        for line in fileinput.FileInput(file,inplace=True,mode='rb'):
                    
            line = line.decode('utf-8')
            #Find the place to put the HTML code for the version dropdown menu
            if '<div role="search">' in line:
                line = line.replace(line,'</div>\n\n\n' + line)
                line = line.replace(line,'    </div>\n' + line)
                line = line.replace(line,'    <div class="select_arrow">\n' + line)
                line = line.replace(line,'    </select>\n' + line)
                #Make the current version the default selection
                for version in tags:
                            
                    href = os.path.relpath(os.path.join('multidocs',version,relpath), os.path.join('multidocs',tag,relpath))
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
            if '<title>DeerLab-development Documentation' in line and tag != 'develop':
                    line = line.replace(line,'<title>DeerLab Documentation</title>')
            if '<li>DeerLab-development Documentation' in line and tag != 'develop':
                    line = line.replace(line,'<li>DeerLab '+tag+' Documentation</li>')    
            if '<h1>DeerLab-development Documentation' in line and tag != 'develop':
                    line = line.replace(line,'<h1>DeerLab '+tag+' Documentation</h1>')  
                    
            #Print line into file
            sys.stdout.buffer.write(line.encode('utf-8'))
            

    #Copy current build into collective docs
    if tag == mostrecent:
        shutil.copytree(docspath, os.path.join(multidocspath,tag))
    else:
        shutil.copytree(docspath, os.path.join(multidocspath,tag))
    
    #Copy the homepage html file to the multidocs build
    if tag == 'develop':
        shutil.copyfile(os.path.join('..','webpage','index.html'), os.path.join('..','multidocs','index.html'))
        shutil.copytree(os.path.join('..','webpage','_static'), os.path.join('..','multidocs','_static'))
        #Go through all lines of index.html
        file = os.path.join('..','multidocs','index.html');
        for line in fileinput.FileInput(file,inplace=True,mode='rb'):
            line = line.decode('utf-8')
            
            #Update the release tag for the HTML links
            if '$RELEASETAG$' in line:
                line = line.replace('$RELEASETAG$',mostrecent)
            #Print into file
            sys.stdout.buffer.write(line.encode('utf-8'))

    print('\n \n VERSION ' + tag + ' COMPLETED                                \n \n', end='\r')

    print('\n')



#Return to the branch the script started at
subprocess.run (["git", "checkout", startbranch])


