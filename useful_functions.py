#useful functions
import numpy as np
import scanpy as sc 
from scipy import stats
import pandas as pd
from progressbar import ProgressBar
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from sklearn.decomposition import NMF


def tfidf_markers(adata, groupby, groupvalue):
    #do TFIDF MARKERS
    dat = adata.X #make the matrix
    target = np.array(adata.obs[groupby].isin([groupvalue])) #boolean of target cluster
    target_len = sum(target) #size of target cluster
    target_dat = dat[target, ] >0
    nobs = np.array(target_dat.sum(axis=0).T) #number of counts for each gene in the cluster
    total_dat = dat > 0
    ntot = np.array(total_dat.sum(axis=0).T) #number of counts for each gene in total
    tf = nobs/sum(target) #term frequency
    universe_len = dat.get_shape()[0] #get the number of cells in total
    idf = np.log(universe_len/ntot) #get inverse document frequency
    score = tf*idf
    q = nobs-1
    m = ntot + universe_len-ntot 
    n = universe_len - ntot
    N = target_len
    pvals = stats.hypergeom.cdf(q, m, n, N) #not sure we are getting this bit totally right.. 
    def FDR(x):
        """
        Assumes a list or numpy array x which contains p-values for multiple tests
        Copied from p.adjust function from R  
        """
        o = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
        ro = [i[0] for i in sorted(enumerate(o), key=lambda v:v[1])]
        q = sum([1.0/i for i in range(1,len(x)+1)])
        l = [q*len(x)/i*x[j] for i,j in zip(reversed(range(1,len(x)+1)),o)]
        l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
        return l
    adj_pvals = np.vstack(FDR(pvals))
    expidf= np.exp(-idf)
    gex = np.array(dat[target, ].mean(axis=0).T) #calculate mean gene expression in the cluster
    gex_out = np.array(dat[~target, ].mean(axis=0).T) #calculate mean gene expression outside the cluster
    gex_glob = np.array(dat.mean(axis = 0).T) #calculate mean gene expression globally
    df = pd.DataFrame(np.column_stack((tf,expidf, gex, gex_out, gex_glob, idf, score, pvals, adj_pvals)), 
            columns = ['GeneFrequency', 'GeneFrequencyGlobal', 'GeneExpressionInCluster','GeneExpressionOutsideCluster',
                       'GeneExpressionGlobal','idf', 'tfidf', 'pval', 'adj_pvals'],
                index = adata.var['Symbol'])
    df = df.sort_values(by = ['tfidf'], ascending = False)
    return df
    
def get_markers(adata, clusters):
    pbar = ProgressBar()
    markers_list = [] #iterate this over a list...
    index_use = list(adata.obs[clusters].cat.categories)
    for i in pbar(list(index_use)):
        mark = tfidf_markers(adata[:,adata.var_names.values[adata.var.highly_variable]], clusters, i)
        markers_list.append(mark)
    cluster_cats = list(adata.obs[clusters].cat.categories)
    markers = dict(zip(cluster_cats, markers_list))
    return markers

def densityplot(x, y):
    # Calculate the point density
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    fig, ax = plt.subplots()
    ax.scatter(x, y, c=z, s=50)
    plt.show()

def run_NMF(adata, n_components):
    model = NMF(n_components=n_components, init='random', random_state=0)
    #sc.pp.highly_variable_genes(adata)
    #mat = adata.X[:, adata.var.highly_variable]
    mat = adata.X
    W = model.fit_transform(mat)
    H = model.components_
    nmf_embedding = pd.DataFrame(W, index = adata.obs_names)
    adata.obsm['X_nmf'] = nmf_embedding
    for i, z in enumerate(W.T):
        adata.obs[f'Z_{i}'] = z
    nmf_loadings = pd.DataFrame(H.T, index = adata.var.Symbol[adata.var.highly_variable])
    return nmf_loadings, adata

def print_loadings(loadings, factor, n_show, adata):
    print('Top loadings by magnitude\n---------------------------------------------------------------------------------------')
    srt = loadings[factor].sort_values(ascending = False)
    print(factor)
    print(srt.head(n_show))
    zs = [f'Z_{factor}']
    sc.pl.umap(adata, color=zs)