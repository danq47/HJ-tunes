# Program to take our top file, and write it in the libtunes format
import sys
from decimal import Decimal
import csv

def convert_to_float(str_num): # top files give numbers as strings with a D instead of E for exponential
    tmp=str_num[:-4]+'E'+str_num[-3:]
    return float(tmp)

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False  

file = sys.argv[1] # takes as input the top file (after rescaling to make it doubly differential) form the command line

f   = open(file, "r")
out = open(file+"-libtunes","w")
out.write("# pT, rapidity, xsection, MC error\n")
out.write("x0\tx1\ttarget\terror\n")

while True:
    x = f.readline() # Loop over lines in file
    if not x: break
    
    line=x.split()

    if line != []:
        line_to_output = "" # declare the line that we will output at the end
        y_mid = 0.0 # declaring the middle of the y bins
        if line[-2] == 'index': # this is a title line
            hist_name = line[1]
            if "yH" in hist_name and "inf" not in hist_name and "ptj1" in hist_name: # Then it is a rapidity slice
                hist_name = hist_name.split("-")
                first = True # this is to get the lower and upper bin edges
                for i in hist_name: # this for loop finds the y slice midpoints
                    if is_float(i):
                        if first:
                            bin_low = float(i)
                            if "minus" in hist_name:
                                bin_low = -1.0 * bin_low
                            first = False
                        else:
                            bin_high = float(i)
                            if "minus" in hist_name:
                                bin_high = -1.0 * bin_high
                y_mid = (bin_high + bin_low)/2.0
# Next, we need to find the pt midpoints for each histogram
        elif "yH" in hist_name and "inf" not in hist_name and "ptj1" in hist_name:
            pt_bin_low = float(line[0])
            pt_bin_high = float(line[1])
            xsec = float(line[2])
            err = float(line[3])
            pt_mid = (pt_bin_high + pt_bin_low)/2.0
            line_to_output = [str(pt_mid),str(y_mid),xsec,err]
            for i in range(2,4):
                line_to_output[i] = '%.8E' % Decimal(str(line_to_output[i]))

            line_to_output="\t".join(line_to_output)
            out.write(line_to_output)
            out.write("\n")

