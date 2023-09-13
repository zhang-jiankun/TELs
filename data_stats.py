import numpy as np
from scipy import stats 

def cdf_function(vector, numbins=25):
    res = stats.cumfreq(vector, numbins=numbins)
    x = res.lowerlimit + np.linspace(0, res.binsize*res.cumcount.size, res.cumcount.size)
    y = res.cumcount/len(vector)
    return x, y


