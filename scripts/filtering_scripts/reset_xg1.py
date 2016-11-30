# reset the bad files in xg1 filtering to check what it looked like before
import glob
import os

files = glob.glob('*-bad-xg1*')
for filename in files:
   filename_tmp=filename.split('-bad-xg1-')
   os.rename(filename,filename_tmp[0])