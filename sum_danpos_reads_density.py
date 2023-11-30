#!/usr/bin/env python

import sys
import numpy as np
import pandas as pd 

matrix = pd.read_csv(sys.argv[1], comment='@', sep='\t', header=None)
header_line = sys.argv[2]
matrix = matrix[range(6, matrix.shape[1])]
vec = matrix.sum().values
np.savetxt(sys.stdout, vec, fmt='%.8f', header=header_line)
