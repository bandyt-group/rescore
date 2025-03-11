import sys
sys.path.append('/home/bmanookian/rescore/')
import rescore_BN as rBN

## User Inputs ##

# If you have a set of trajectories. Provide the names and the lengths of the trajs
files=None    #['Rep1.csv','Rep2.csv','Rep3.csv']
lengths=None   #[100,100,100]

# If you have a single file that has already combned the trajectroies use filename
filename='filename.csv'

# Provide the .dot file for your BN universal graph
dotfile='./graph.dot'

# Data 'discrete' or 'continuous'
datatype='continuous'

# Output file for rescore values
outfile='out.csv'

# nuber of cores
ncore=20
##


## Run codes ##
# Use the folloiwng for sets of trajectories
Rescore=rBN.BN_Rescore(dotfile=dotfile,files=files,data=filename,lengths=lengths)
Rescore.rescore(datatype=datatype,ncore=ncore)
Rescore.table_write(outfile)

##
