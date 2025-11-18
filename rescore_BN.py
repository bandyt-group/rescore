import numpy as np
from scipy.spatial import cKDTree
from scipy.special import digamma
import scipy.spatial as ss
import networkx as nx
from networkx.drawing.nx_pydot import read_dot
from networkx.algorithms.moral import moral_graph
from functools import partial
import os
import csv
import itertools
import multiprocessing
#from entropy_estimators_continuos import condMutInf_mixed_multid

## Data Loader ##
def csvreader(inputfile):
    out = []
    with open(inputfile, newline = '') as file:
        reader = csv.reader(file)
        for i,row in enumerate(reader):
            out.append(row)
    return np.array(out)

def csvwriter(outfile,data):
    with open(outfile, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        for row in data:
            csv_writer.writerow(row)

## Entropy and Mutual Information Functions **See bottom of code for Cont MI**  ##

def encode(c):
    try:
        b=np.ones(c.shape[1],dtype=int)
    except Exception:
        c=np.column_stack(c)
        b=np.ones(c.shape[1],dtype=int)
    b[:-1]=np.cumprod((c[:,1:].max(0)+1)[::-1])[::-1]
    return np.sum(b*c,1)


def mi(A):
    X,Y=A
    return H(X)+H(Y)-joinH((X,Y))

def H(i):
        """entropy of labels"""
        p=np.unique(i,return_counts=True)[1]/i.size
        return -np.sum(p*np.log2(p))

def joinH(i):
    pair=np.column_stack((i))
    en=encode(pair)
    p=np.unique(en,return_counts=True)[1]/len(en)
    return -np.sum(p*np.log2(p))

## Parallel Function ##

def runParallel(foo,iter,ncore):
    pool=multiprocessing.Pool(processes=ncore)
    try:
        out=(pool.map_async( foo,iter )).get()
    except KeyboardInterrupt:
        print ("Caught KeyboardInterrupt, terminating workers")
        pool.terminate()
        pool.join()
    else:
        #print ("Quitting normally core used ",ncore)
        pool.close()
        pool.join()
    try:
        return out
    except Exception:
        return out

def readDot(f):
    return nx.DiGraph(read_dot(f))

def knn(data, k):
    tree = cKDTree(data)
    distances, indices = tree.query(data, k=k+1)
    return distances[:,1:], indices[:,1:]

def getlabelarray(lengths):
    k=np.array([np.sum(lengths[:i]) for i in range(len(lengths)+1)]).astype(int)
    int_indx=[(k[i],k[i+1]) for i in range(len(k)-1)]
    q=np.ones(k[-1])
    for i in range(len(int_indx)):
        q[int_indx[i][0]:int_indx[i][1]]=i
    return q

def make_batches(total_len, batch_size):
    return [(i, min(i + batch_size, total_len)) for i in range(0, total_len, batch_size)]

def getunion(outs):
    U=np.sort(np.unique(np.concatenate([out[0] for out in outs])))
    return [aligndata(out,U) for out in outs]
    
def aligndata(data,un_labs):
    m=np.zeros([data[1:].shape[0],un_labs.shape[0]]).T
    for i,k in enumerate(np.where(np.in1d(un_labs,data[0])==True)[0]):
        m[k]=data[1:][:,i]
    return np.vstack((un_labs,m.T))

def combinecsvfiles(files,union=False):
    outs=[csvreader(f) for f in files]
    if union:
        outs=getunion(outs)     
    concat=[]
    lengths=[]
    for i,out in enumerate(outs):
        lengths.append(out[1:].shape[0])
        if i==0:
            concat=out[:,np.argsort(out[0])]
        if i!=0:
            concat=np.vstack((concat,out[1:,np.argsort(out[0])].astype(float)))
    indx=np.split(np.arange(0,concat[1:].shape[0]),np.cumsum(lengths[:-1]))
    return concat,indx

class BN_Rescore:
    def __init__(self,dotfile,data=None,files=None,union=False,subsets=None,discrete=False):
        ## combine csv files ##
        if files is not None:
            data,self.subset_indx=combinecsvfiles(files,union=union)
        ## Load file
        if type(data)==str:
            data=csvreader(data)

        ## discrete or continuous
        self.discrete=discrete

        ## Subset labels ##
        #self.subset_indx=subsets
        if data[0,-1]=='label':
            self.subset_labels=data[1:,-1]
            self.subset_label_set=np.sort(list(set(self.subset_labels)))
            self.subset_indx=[np.where(self.subset_labels==i)[0] for i in self.subset_label_set] 
        if type(subsets)==list:
            self.subset_indx=subsets

        ## Input variables and data with sorting ##
        self._input_data=data[1:]
        if discrete is not True:
            self._input_data=self._input_data.astype(float)
        if discrete:
            self._input_data=self._input_data.astype(float).astype(int)
        self._input_variables=data[0]
        i_sort=np.argsort(self._input_variables)
        self.variables=self._input_variables[i_sort]
        self.data=self._input_data[:,i_sort]
        self.n_cells=self.data.shape[0]
        
        ## Graph input from dotfile and reading nodes and edges ##
        self.G=readDot(dotfile)
        self.nodes=np.array(list(self.G.nodes()))
        self.edges=np.array(list(self.G.edges()))
        self.varsNotInG=np.setdiff1d(self.variables,self.nodes)

        if len(self.varsNotInG)>0:
            i=~np.in1d(self.variables,self.varsNotInG)
            self.data=self.data[:,i]
            self.variables=self.variables[i]
        self.mapVarNodes=np.searchsorted(self.variables,self.nodes)
        self.n_nodes=len(self.G.nodes)
        self.n_variables=len(self.variables)
        self.G_moral=moral_graph(self.G)
        self.map_variable=dict(zip(self.variables,np.arange(self.variables.size)))

## RESCORING function ##
    def rescore(self,ncore=1):
        if self.discrete:
            self.rescored_out=self.disc_MI_on_subsets(ncore=ncore)
            return
        self.rescored_out=self.cont_MI_on_subsets(ncore=ncore)
        return


## Discrete data MI functions ##
    def mi_on_edge(self,u,v):
        if u not in self.map_variable:
            return 0
        if v not in self.map_variable:
            return 0
        x=self.data[:,self.map_variable[u]]
        y=self.data[:,self.map_variable[v]]
        return mi([x,y])
    def _mi_on_edge_subset_data(self,edge,subset):
        u,v=edge
        x=self.data[subset,self.map_variable[u]]
        y=self.data[subset,self.map_variable[v]]
        return mi([x,y])
    def disc_MI_on_edges(self,subset,edges):
        return np.array([self._mi_on_edge_subset_data(edge=ed,subset=subset) for ed in edges])
    def disc_MI_on_batch(self,batch,subsets,edges):
        return [self.disc_MI_on_edges(subset=s,edges=edges) for s in subsets[batch[0]:batch[1]]]
    def disc_MI_on_subsets(self,edges=None,subsets=None,moral=False,batches=None,ncore=4):
        if edges is None:
            edges=self.G.edges
        if moral:
            edges=self.G_moral.edges
        if self.subset_indx is None:
            subset=np.arange(self.n_cells)
            func=partial(self._mi_on_edge_subset_data,subset=subset)
            return runParallel(func,edges,ncore)
#            return self.disc_MI_on_edges(edges=edges,subset=subset,ncore=ncore)
#        if subsets == 'all':
        subsets=self.subset_indx
        if batches is not None:
          func=partial(self.disc_MI_on_batch,edge_list=edges)
          return runParallel(func,batches,ncore=ncore)
        func=partial(self.disc_MI_on_edges,edges=edges)
        return runParallel(func,subsets,ncore)
#        return [self.disc_MI_on_edges(edges=edges,subset=s,ncore=ncore) for s in subsets]
    #def disc_MI_on_batch(self,subsets,edges,batchs):
    def disc_MI_on_batch(self, batch, edge_list):
        return [self.disc_MI_on_edges(subset=s,edges=edge_list) for s in self.subset_indx[batch[0]:batch[1]]]
### Continuous MI Functions ###

    def _mi_Mixed_KSGm_on_edge(self,edge,k_bin=5):
        u,v=edge
        x=self.data[:,self.map_variable[u]]
        y=self.data[:,self.map_variable[v]]
        if len(set(x))==1:
            return 0.0
        if len(set(y))==1:
            return 0.0
        return Mixed_KSGm(x,y,k_bin)
    def _mi_Mixed_KSGm_on_edge_subset_data(self,edge,subset,k_bin=5):
        u,v=edge
        x=self.data[subset,self.map_variable[u]]
        y=self.data[subset,self.map_variable[v]]
        if len(set(x))==1:
            return 0.0
        if len(set(y))==1:
            return 0.0
        return Mixed_KSGm(x,y,k_bin)
    def cont_MI_on_edges(self,subset,edges,k_bin=5,ncore=1):
        return np.array([self._mi_Mixed_KSGm_on_edge_subset_data(edge=ed,subset=subset) for ed in edges])
#        func=partial(self._mi_Mixed_KSGm_on_edge_subset_data,subset=subset,k_bin=k_bin)
#        return runParallel(func,edges,ncore)
    def cont_MI_on_subsets(self,edges=None,moral=False,k_bin=5,ncore=1):
        if edges is None:
            edges=self.G.edges
        if moral:
            edges=self.G_moral.edges
        if self.subset_indx is None:
            subset=np.arange(self.n_cells)
            func=partial(self._mi_Mixed_KSGm_on_edge_subset_data,subset=subset,k_bin=k_bin)
            return runParallel(func,edges,ncore)
#            return self.cont_MI_on_edges(edges=edges,subset=subset,ncore=ncore)
        subsets=self.subset_indx
        func=partial(self.cont_MI_on_edges,edges=edges,ncore=ncore)
        return runParallel(func,subsets,ncore)
#        return [self.cont_MI_on_edges(edges=edges,subset=s,ncore=ncore) for s in subsets]
    
## Spatial Continuous Data Functions ##
    def _mi_Mixed_KSGm_spatial(self,kernal,edge,k_bin=5):
        u,v=edge
        x=self.data[kernal,self.map_variable[u]]
        y=self.data[kernal,self.map_variable[v]]
        if len(set(x))==1:
            return 0.0
        if len(set(y))==1:
            return 0.0
        return Mixed_KSGm(x,y,k_bin)

    def run_kernals(self,edge,kernals):
        return np.array([self._mi_Mixed_KSGm_spatial(k,edge) for k in kernals])

    def cont_MI_on_subset(self,edge,kernals,k_bin=5,ncore=20):
        func=partial(self._mi_Mixed_KSGm_spatial,edge=edge)
        return runParallel(func,kernals,ncore=ncore)

    def cont_MI_on_spatial(self,edges=None,kernals=None,moral=False,k_bin=5,ncore=1):
        if edges is None:
            edges=self.edges
        if moral:
            edges=self.G_moral.edges
        func=partial(self.run_kernals,kernals=kernals)
        return np.array(runParallel(func,edges,ncore))

## Writing Rescored to output table ##
    def table(self):
        self.table=self.edges
        if self.subset_indx is None:
            self.table=np.column_stack((self.table,self.rescored_out))
            return
        for i,s in enumerate(self.subset_indx):
            self.table=np.column_stack((self.table,self.rescored_out[i]))
    def table_write(self,outfile):
        self.table()
        if self.subset_indx is None:
            subset_labels=['rescored']
        if self.subset_indx is not None: 
            subset_labels=np.array([f'subset{i+1}' for i in range(len(self.subset_indx))])
        title_row=np.insert(subset_labels,0,np.array(['source','target']))
        csvwriter(outfile,np.vstack((title_row,self.table)))

## Weighted degree ##
    def _get_map2adj(self,is_moral=False):
            edges=self.G.edges
            if is_moral:
                edges=self.G_moral.edges
            return tuple(zip(*[[self.map_variable[u],self.map_variable[v]] for u,v in  edges]))
    def conpute_wd(self,edge_scores,map2adj):
        A=np.zeros((self.n_variables,self.n_variables))
        A[map2adj]=edge_scores
        return A.sum(1)

    def conputeAll_wd(self, edge_scores=None,is_moral=False, foo=None,**args):
        """"" return wd for each subset""""" 
        if edge_scores is None:
            print('computing edges_scores')
            edge_scores=foo(**args)
        map2adj=self._get_map2adj(is_moral)
        return np.array([self.conpute_wd(scores,map2adj) for scores in edge_scores])

## Partial MI functions ##
    def _pmi_on_node_all_subsets(self,node,k_bin,subsets=None):
        all_pa=list(self.G.predecessors(node))
        n_pa=len(all_pa)
        if len(all_pa)==0:
            return None
        data_pa=np.array([self.data[:,self.map_variable[u]] for u in all_pa if u in self.map_variable]).T
        data_child=self.data[:,self.map_variable[node]]
        if subsets is None:
            subsets=[np.arange(self.n_cells)]
        if n_pa == 1:
            pmi=[Mixed_KSGm(data_child[subset],data_pa[subset],k_bin) for subset in subsets]
            return [[(all_pa[0],node),pmi]]
        if n_pa == 2:
            pmi=np.array([[condMutInf_mixed_multid(k_bin,data_child[subset],data_pa[subset,i],self.data[subset,j]) for i,j in ((0,1),(1,0)) ] for subset in subsets]).T
            return [[(pai,node),pmi_i] for pai,pmi_i in zip(all_pa,pmi)]
        all_parents_but_one=list(itertools.combinations(data_pa.T,n_pa-1))
        pmi=[[condMutInf_mixed_multid(k_bin,data_child[subset],data_pa[subset,i],np.column_stack(Z)[subset]) for i,Z in enumerate(all_parents_but_one) ] for subset in subsets]

        return [[(pai,node),pmi_i] for pai,pmi_i in zip(all_pa[::-1],pmi)]

    def pmi_bn(self,k_bin=3,subsets=None,ncore=-1):
        foo_pmi=partial(self._pmi_on_node_all_subsets,k_bin=k_bin,subsets=subsets)
        if ncore==-1:
            ncore=0
            if subsets is not None:
                ncore=np.min([len(subsets),os.cpu_count()])
        if ncore==0:
           pmi = [foo_pmi(node) for node in self.G.nodes]
        else:
           pmi = runParallel(foo_pmi,self.G.nodes,ncore)
        tmp_dict={}
        for x in pmi:
            if x is not None:
                for a,b in x:
                    tmp_dict.update({a:b})
       
        pmi=np.array([tmp_dict[edge]  for edge in self.G.edges])
        if pmi.shape[1]==1:
            pmi=pmi.flatten()
        return pmi
    def update_graph(self,G=None,nodes_property=None,edges_property=None,node_property_name='prop',edge_property_name='score'):
        if G is None:
            G=self.G
        if G =='moral':
            G=self.G_moral
        update_graph(G,nodes_property=nodes_property,edges_property=edges_property,node_property_name=node_property_name,edge_property_name=edge_property_name)
        
def update_graph(G,nodes_property=None,edges_property=None,node_property_name='prop',edge_property_name='score'):
    if nodes_property is not None:
        update_nodes_graph(G,nodes_property,node_property_name)
    if edges_property is not None:
        update_edges_graph(G,edges_property,edge_property_name)
def update_nodes_graph(G,nodes_property,node_property_name='prop'):
    if type(nodes_property) != dict:
        nodes_property=dict(zip(G.nodes,nodes_property))
    nx.set_node_attributes(G, nodes_property, name=node_property_name)
def update_edges_graph(G,edges_property,edge_property_name='score'):
    if type(edges_property) != dict:
        edges_property=dict(zip(G.edges,edges_property))
    nx.set_edge_attributes(G, edges_property, name=edge_property_name)
    
   
def Mixed_KSGm(x,y,k,onlyPos=True):
     '''
             Estimate the mutual information I(X;Y) of X and Y from samples {x_i, y_i}_{i=1}^N
             Using *Mixed-KSG* mutual information estimator
 
             Input: x: 2D array of size N*d_x (or 1D list of size N if d_x = 1)
             y: 2D array of size N*d_y (or 1D list of size N if d_y = 1)
             k: k-nearest neighbor parameter
 
             Output: one number of I(X;Y)
     '''
     assert len(x)==len(y), "Lists should have same length"
     assert k <= len(x)-1, "Set k smaller than num. samples - 1"

     try:
         N = len(x)
         if x.ndim == 1:
                 x = x.reshape((N,1))
         dx = len(x[0])
         if y.ndim == 1:
                 y = y.reshape((N,1))
         dy = len(y[0])
         data = np.concatenate((x,y),axis=1)
     
         tree_xy = ss.cKDTree(data)
         tree_x = ss.cKDTree(x)
         tree_y = ss.cKDTree(y)
     
         knn_dis=np.array(tree_xy.query(data,k+1,p=float('inf'))[0][:,k])
         j=np.where(knn_dis==0)[0]
         if len(j)==0:
             knn_dis-=1e-15
             nx,ny=zip(*[[len(tree_x.query_ball_point(x[i],knn_dis[i],p=np.inf)),len(tree_y.query_ball_point(y[i],knn_dis[i],p=np.inf))] for i in range(N)])
             I=digamma(k)+np.log(N)-np.mean(digamma(nx)+digamma(ny))
         else:
             knn_dis-=1e-15
             kh=np.zeros(N)+digamma(k)
             nx,ny=zip(*[[len(tree_x.query_ball_point(x[i],knn_dis[i],p=np.inf)), len(tree_y.query_ball_point(y[i],knn_dis[i],p=np.inf))] for i in range(N) if i not in j])
             ck=[len(z) for z in tree_xy.query_ball_point(data[j],1e-15,p=np.inf)]
             cx=[len(z) for z in tree_x.query_ball_point(x[j],1e-15,p=np.inf)]
             cy=[len(z) for z in tree_y.query_ball_point(y[j],1e-15,p=np.inf)]
             kh[-j.size:]=digamma(ck)
             nx=list(nx)+cx
             ny=list(ny)+cy
             I=np.log(N)+np.mean(kh-digamma(nx)-digamma(ny))
         if (onlyPos)&(I<0):
             I=0.0 
         return I


     except ValueError:
         print('Making discrete')
         not_discrete_x=np.sum(np.unique(x,return_counts=True)[1]==1)
         not_discrete_y=np.sum(np.unique(y,return_counts=True)[1]==1)
         if (not_discrete_x<=1)|(not_discrete_x<=2):
             return mi((x,y))
         print ("check for error\nu\n: ",u,x,'\nv:\n',v,y)


