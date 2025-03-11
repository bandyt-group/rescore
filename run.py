import sys
sys.path.append('/home/bmanookian/rescore/')
import rescore_BN as rBN

## User Inputs ##
files=['Rep1.csv','Rep2.csv','Rep3.csv']
lengths=[100,100,100]
dotfile='./graph.dot'
datatype='continuous'
ncore=20
outfile='out.csv'
##


## Run codes ##
Rescore=rBN.BN_Rescore(dotfile=dotfile,lengths=lengths,files=files)
Rescore.rescore(datatype=datatype,ncore=ncore)
Rescore.table_write(outfile)
##
