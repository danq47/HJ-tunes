# reset the bad files in xg1 filtering to check what it looked like before
import glob
import os

files = glob.glob('*dat-save')
for filename in files:
	filename_tmp=filename.split('-')
	os.rename(filename,'-'.join(filename_tmp[0:-1]))