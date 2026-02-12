import pandas as pd
import numpy as np
import graphtools as gt
import phate
import scprep
import meld
#import cmocean
import sklearn
import scipy
import seaborn as sns


# making sure plots & clusters are reproducible
np.random.seed(42)

import scanpy as sc
adata = sc.read_h5ad('/cfs/klemming/home/s/sumalw/project/meld_vfc/SSRI.h5ad')

data= adata.X
data = scprep.filter.filter_rare_genes(data)
data_libnorm, libsize = scprep.normalize.library_size_normalize(data, return_library_size=True)
data_sqrt = np.sqrt(data_libnorm)
data_pca = scprep.reduce.pca(data_sqrt) ## take a long time
np.save('/cfs/klemming/home/s/sumalw/project/meld_vfc/sert_data_all_pca.npy', data_pca)
benchmarker = meld.Benchmarker()
benchmarker.fit_phate(data_pca);

from joblib import Parallel, delayed

def simulate_pdf_calculate_likelihood(benchmarker, seed, beta):
    benchmarker.set_seed(seed)
    benchmarker.generate_ground_truth_pdf()

    benchmarker.generate_sample_labels()
    benchmarker.calculate_MELD_likelihood(beta=beta)
    MELD_mse = benchmarker.calculate_mse(benchmarker.expt_likelihood)
    return MELD_mse, seed, beta, benchmarker.graph.knn

knn_range = np.arange(1,25)
beta_range = np.arange(1,200)

results = []

with Parallel(n_jobs=16) as p:
    for knn in knn_range:
        # doing this outside the parallel loop because building the graph takes the longest
        benchmarker.fit_graph(adata.X, knn=knn)
        print(knn)
        curr_results = p(delayed(simulate_pdf_calculate_likelihood)(benchmarker, seed, beta) \
                                       for seed in range(25) for beta in beta_range)
        curr_results = pd.DataFrame(curr_results, columns = ['MSE', 'seed', 'beta', 'knn'])
        #results.append(curr_mse)
        results.append(curr_results)

results = pd.concat(results, axis=0)

results.to_csv('/cfs/klemming/home/s/sumalw/project/meld_vfc/Sert_meld_optimal_results.csv')
