import pandas as pd
import linecache    # Allows us to read a certain line of the file
import glob
import os

seeds=[]
xsec=[]
err=[]
powheg_ave=[]
running_pwg_err=[]
weighted_ave=[]
weighted_err=[]
for filename in glob.iglob('pwg-st2-*-stat.dat'):
	xsec_and_err=linecache.getline(filename, 5) # this is the line containing the error and xsec
	xsec_and_err=xsec_and_err.split() 	# Split on any whitespace 
	xsec.append(float(xsec_and_err[6]))		# The xsec
	err.append(float(xsec_and_err[8]))			# The error
	filename=filename.split("-")
	seeds.append(filename[2])			# The seed (from the filename)	
d = {'Seed' : seeds , 'Xsec' : xsec , 'Error' : err}
df = pd.DataFrame(d, columns=['Seed','Xsec','Error']) # Make a dataframe out of the seeds	
y=0.0
error=0.0
ave=0.0
err_pwg=0.0
for i in range(df['Xsec'].size):		# Work out running averages
	ave = (df['Xsec'][i] + ave*i)/(i+1)
	y = y + df['Xsec'][i]/(df['Error'][i]**2)
	error = error + 1.0/(df['Error'][i]**2)
	weighted_ave.append(y/error)
	powheg_ave.append(ave)
	weighted_err.append(1.0/(error**0.5))	
	err_pwg = err_pwg + df['Error'][i]**2
	running_pwg_err.append((err_pwg/((i+1)**2))**0.5)

df['POWHEG mean']=powheg_ave
df['POWHEG running error']=running_pwg_err
df['Weighted Average']=weighted_ave   
df['Weighted Error']=weighted_err	

print(df.head(80))
