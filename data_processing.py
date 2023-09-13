def read_matrix_score(matrix, norm=False, cnorm=False):
    """ Read `deeptools computeMatrix` output """
    df = pd.read_csv(matrix, sep='\t', comment='@', header=None)
    df = df.drop(columns=list(range(6)))
    
    if norm:
        df[df>np.nanquantile(df, 0.99)] = np.nanquantile(df, 0.99)
    
    return df.mean(axis=0).values

def read_bw_score(filename, sample=None):
    """ Read `bwtool summary` output """

    df = pd.read_csv(filename, sep='\t')
    df = df.rename(columns={'mean': 'signal'})    
    if sample:
        df = df.sample(n=sample, random_state=111)
    df = df.sort_values(['signal'], ascending=False)
    
    return df


def filter_insulation_boundaries(insul_path, window=100000, threshold=0.2):
    """ 
        Filter weak insulation boundaries.
        See https://github.com/open2c/cooltools/issues/64.
    """
    insul_df = pd.read_csv(insul_path, sep='\t')
    insul_df = insul_df[(insul_df['is_bad_bin']==False) & (insul_df[f'is_boundary_{window}']==True)]
    insul_df = insul_df[insul_df[f'boundary_strength_{window}']>threshold] 
    insul_df = insul_df[['chrom', 'start', 'end']]
    return inful_df

