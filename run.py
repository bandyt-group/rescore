import sys
sys.path.append('/home/bmanookian/rescore/')
import rescore_BN as rBN

## User Inputs ##
# Number of cores #
ncore=20

## Graph input in the form of a dotfile ##
dotfile='/path/to/dotfile'


## File/Data input ##

# If you have a set of trajectories. Provide the names and the lengths of the trajs
#    otherwise put None
files=None    #['Rep1.csv','Rep2.csv','Rep3.csv']

# If you have a single file that has already combned the trajectroies use filename
filename='filename.csv'

# Data 'discrete' or 'continuous'
datatype='continuous'

## Output ##
# Output file for rescore values
outfile='out.csv'

##


## Run codes ##
# Use the folloiwng for sets of trajectories
Rescore=rBN.BN_Rescore(dotfile=dotfile,data=filename,files=files)
Rescore.rescore(datatype=datatype,ncore=ncore)
Rescore.table_write(outfile)

##
