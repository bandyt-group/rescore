import numpy as np
from scipy.special import digamma
import scipy.spatial as ss
import igraph as ig


def entropy(i):
    """entropy of labels"""
    p=np.unique(i,return_counts=True)[1]/i.size
    return -np.sum(p*np.log(p))
def mutualInformation(X,Y):
    return entropy(X)+entropy(Y)-jointEntropy((X,Y))

def jointEntropy(i):   
   pair=np.column_stack((i))
   en=encode(pair)
   p=np.unique(en,return_counts=True)[1]/len(en)
   return -np.sum(p*np.log(p))
def encode(c):
     b=np.ones(c.shape[1],dtype=int)
     b[:-1]=np.cumprod((c[:,1:].max(0)+1)[::-1])[::-1]
     return np.sum(b*c,1)



def Mixed_KSG(input):
    '''
		Estimate the mutual information I(X;Y) of X and Y from samples {x_i, y_i}_{i=1}^N
		Using *Mixed-KSG* mutual information estimator

		Input: x: 2D array of size N*d_x (or 1D list of size N if d_x = 1)
		y: 2D array of size N*d_y (or 1D list of size N if d_y = 1)
		k: k-nearest neighbor parameter

		Output: one number of I(X;Y)
    '''

    x,y,k=input
    assert len(x)==len(y), "Lists should have same length"
    assert k <= len(x)-1, "Set k smaller than num. samples - 1"
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

    knn_dis = [tree_xy.query(point,k+1,p=float('inf'))[0][k] for point in data]
    ans = 0
    for i in range(N):
        kp, nx, ny = k, k, k
        if knn_dis[i] == 0:
            kp = len(tree_xy.query_ball_point(data[i],1e-15,p=float('inf')))
            nx = len(tree_x.query_ball_point(x[i],1e-15,p=float('inf')))
            ny = len(tree_y.query_ball_point(y[i],1e-15,p=float('inf')))
        else:
            nx = len(tree_x.query_ball_point(x[i],knn_dis[i]-1e-15,p=float('inf')))
            ny = len(tree_y.query_ball_point(y[i],knn_dis[i]-1e-15,p=float('inf')))
        ans += (digamma(kp) + np.log(N) - digamma(nx) - digamma(ny))/N
    return ans


def get_sub_space_indx(X):
    sub_space_indx=[]
    indx=0
    for x in X:
        sub_indx=[indx]
        if x.ndim>1:
            sub_indx=[indx+j for j in range(x.shape[1])]
        indx=sub_indx[-1]+1
        sub_space_indx.append(sub_indx)
    return sub_space_indx


def GDM(X,G,k):
    N = len(X)
    tree_X = ss.cKDTree(X)
    knn_dis = np.array([tree_X.query(point,k+1,p=float('inf'))[0][k] for point in X])
    id0=knn_dis==0
    digamma_k=digamma(k)
    
    if np.any(id0):
        k_tilde = np.array([len(tree_X.query_ball_point(x,1e-15,p=float('inf'))) for x in X[id0]])
        digamma_k=np.repeat(digamma_k,N)
        digamma_k[id0]=digamma(k_tilde)
    inq,z=zip(*[inquire_n(X,k,node,G.predecessors(node),knn_dis) for node  in range(G.vcount())])
    return np.mean(digamma_k+np.sum(inq,0))+(np.sum(z)-1)*np.log(N)

def GDM_multidim(X,G,k):
    sub_space_indx=get_sub_space_indx(X)
    X=np.column_stack(X)
    N = len(X)
    tree_X = ss.cKDTree(X)
    knn_dis = np.array([tree_X.query(point,k+1,p=float('inf'))[0][k] for point in X])
    id0=knn_dis==0
    digamma_k=digamma(k)
    if np.any(id0):
        k_tilde = np.array([len(tree_X.query_ball_point(x,1e-15,p=float('inf'))) for x in X[id0]])
        digamma_k=np.repeat(digamma_k,N)
        digamma_k[id0]=digamma(k_tilde)
    #return G,knn_dis,sub_space_indx,X,k 
    inq,z=zip(*[inquire_n(X,k,sub_space_indx[node],G.predecessors(node),knn_dis) for node  in range(G.vcount())])
    return np.mean(digamma_k+np.sum(inq,0))+(np.sum(z)-1)*np.log(N)

def inquire_n(X,k,l,pa,knn_dis):
    def merge_dim(x,y):
        if np.ndim(x)==0:
            x=[x]
        return x+y
    if pa:
        return digamma(count_nbh_subset(X[:,pa],knn_dis,k))-digamma(count_nbh_subset(X[:,merge_dim(l,pa)],knn_dis,k)),0
        #return np.log(count_nbh_subset(X[:,pa],knn_dis,k)+1)-np.log(count_nbh_subset(X[:,[l]+pa],knn_dis,k)+1),0
    return -digamma(count_nbh_subset(X[:,l],knn_dis,k)),1
    #return -np.log(count_nbh_subset(X[:,l],knn_dis,k)+1),1
    
def count_nbh_subset(S,knn_dis,k):
    if S.ndim==1:
        S=S.reshape(S.size,1)
    tree_n_pa = ss.cKDTree(S)
    n_subset = np.array([count_nbh_i(tree_n_pa,s,d) for s,d in zip(S,knn_dis)])
    return n_subset
def count_nbh_i(tree_n_pa,s,d):
    if d==0:
        return  len(tree_n_pa.query_ball_point(s,1e-15,p=float('inf')))
    return len(tree_n_pa.query_ball_point(s,d-1e-15,p=float('inf')))

    
def total_cor_mixed(k,data):
    g=ig.Graph(directed=True)
    g.add_vertices(data.shape[1])
    return GDM(data,g,k)
def condMutInf_mixed(k,data):
    """""in data first column is x remaining are child of x"""""
    g=ig.Graph(directed=True)
    g.add_vertices(data.shape[1])
    [g.add_edge(0,i) for i in range(1,data.shape[1])]
    return GDM(data,g,k)
def condMutInf_mixed_multid(k,x,y,z):
    """""not working PMI(x;y|z)"""""
    g=ig.Graph(directed=True)
    g.add_vertices(3)
    g.add_edges([(2,0),(2,1)])
# GDM need to handle not homogeneous dimension
    return GDM_multidim((x,y,z),g,k)


def Mixed_KSGm_safe(x,y,k,onlyPos=True):
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
         not_discrete_x=np.sum(np.unique(x,return_counts=True)[1]==1)
         not_discrete_y=np.sum(np.unique(y,return_counts=True)[1]==1)
         if (not_discrete_x<=1)|(not_discrete_x<=2):
             return mutualInformation(x,y)
         print ("check for error\nu\n: ",u,x,'\nv:\n',v,y)

