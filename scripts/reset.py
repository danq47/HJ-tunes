import os
import glob

for filename in glob.iglob('*.dat-save'):
	filename=filename.split('-')
	filename=filename[:4]
	filename='-'.join(filename)
	os.rename(filename+'-save',filename)