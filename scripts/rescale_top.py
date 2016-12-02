import sys

def convert_to_float(str_num): # top files give numbers as strings with a D instead of E for exponential
    tmp=str_num[:-4]+'E'+str_num[-3:]
    return float(tmp)

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False    

file = sys.argv[1]          # takes as input a top file

f   = open(file, "r")
out = open(file+"-rescaled","w")

bin_high=0
bin_low=0

while True:
    x = f.readline() # Loop over lines in file
    if not x: break
    
    line=x.split()
    
    if line != []:
        if line[-2] == 'index': # find the title line for each hist
            hist_name=line[1]
            if "yH" in hist_name and "inf" not in hist_name: # then it is one that is sliced in rapidity
                hist_name=hist_name.split("-")
                first=True # to get the two different bin edges
                for i in hist_name:
                    if is_float(i):
                        if first:
                            bin_low = float(i)
                            first = False
                        else:
                            bin_high = float(i)
                slice_width = abs( bin_low - bin_high )
            else:
                slice_width = 1.0
        else: # so it's not a title line and it's not a blank line, therefore we divide column 3 and 4 by slice width
            for i in range(2,4):
                line[i] = str(convert_to_float(line[i])/slice_width)

        line = " ".join(line) # remake the line in the .top format for writing
    else:
        line = ""
        
    out.write(line)
    out.write("\n")

                