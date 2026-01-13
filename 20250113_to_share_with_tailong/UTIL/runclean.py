#!/home/disk/eos1/mwyant/anaconda3a/bin/python 
"""
Created on Tue Nov 17 14:22:45 2015

@author: Matt Wyant
Delete all SAM output files for a run
runclean.py [casename] [prm file name]

To be run from the main SAM directory. Will offer to delete files from the OUT_* and
and RESTART directories. 

The ./[casename]/[prm file name] file is opened to locate the caseid. The caseide is 
used to construct the potential output/restart filenames.
  
The second argument is optional, without it the script will try to open
the ./[casename]/prm file.

"""

import sys
import os
import glob
import re

nargs = len(sys.argv)

if nargs < 2 or nargs > 3:
    raise Exception('Incorrect number of arguments, ' + \
    'should be runclean.py [casename] [prm file name]')
    
casename = sys.argv[1]    
        
if nargs == 2:
    prmfilename = 'prm'
else:
    prmfilename = sys.argv[2]

prmpath = '%s/%s' % (casename, prmfilename)
print('\nRunning runclean......')
print('casename: ' + casename)  
print('prm filename: ' + prmfilename)

# open the file, find the caseid from inside the prm file   
f = open(prmpath, 'r')
caseid = ''
filetext = f.read()
for line in filetext.splitlines():  
    ind_quote = line.find("'", line.find('=', line.find('caseid')))
    if ind_quote !=1:       
        ind_endquote = line.find("'", ind_quote+1)
        if ind_endquote != -1:
            caseid = line[ind_quote+1:ind_endquote]
            break
    
f.close()     
print('caseid: ' + caseid + '\n')   
if caseid == '':
    raise Exception('No caseid found in %s' % prmfilename)


# identify all files to delete
outroot = '%s_%s' % (casename, caseid)

candidates = glob.glob('OUT_*/%s*' % outroot)
candidates += glob.glob('RESTART/%s*' % outroot)

outfilepatterns = ['OUT_STAT/%s\.nc', 'OUT_STAT/%s\.stat', \
    'OUT_2D/%s_[0-9]+\.2Dbin', 'OUT_3D/%s_[0-9]+\.bin2D', \
    'OUT_3D/%s_[0-9]+_[0-9]{10}\.bin3D', \
    'RESTART/%s_misc_restart\.bin' , 'RESTART/%s_[0-9]+_restart\.bin',\
    'RESTART/%s_[0-9]+_restart_rad\.bin']
filepattern = '|'.join([s % outroot for s in outfilepatterns])

files_to_delete = [];
for c in candidates:
    if re.match(filepattern, c):
        files_to_delete.append(c)
        

# ask the user
print('Found %s output files:\n' % len(files_to_delete))
for f in files_to_delete: print(f)

# delete the files
delete_str = input('\nDelete them now (y/[n]) ??? ')
if delete_str == 'y':
    print('Deleting %d files...\n' % len(files_to_delete))
    for f in files_to_delete:
        cmdstr = 'rm ' + f
        print(cmdstr)  
        os.system(cmdstr)
else:
    print('No files deleted.')        
        

	


  
