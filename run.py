import sys
sys.path.append('/home/bmanookian/rescore/')
import rescore_BN as rBN

## User Inputs ##
# Number of cores (Hint: Use number of subsets)#
ncore=3

## Graph input in the form of a dotfile ##
dotfile='/path/to/dotfile'


## File/Data input ##

# If you have a set of trajectories. Provide the names in a list, otherwise put None
files=None    #['Rep1.csv','Rep2.csv','Rep3.csv']

# If the labels need to be combined via union
union=False

# If you have a single file that has already combined the trajectroies use filename
# Make sure to have the final column 'label' for subset labels
# ** If no label column need to provide subset indices as input **
filename='filename.csv'
subsets=None

# Datatype 'discrete'(True) or 'continuous'(False)
discrete=True

## Output ##
# Output file for rescore values
outfile='out.csv'

##


## Run codes ##
Rescore=rBN.BN_Rescore(dotfile=dotfile,data=filename,files=files,union=union,subsets=subsets,discrete=discrete)
Rescore.rescore(ncore=ncore)
Rescore.table_write(outfile)

##
