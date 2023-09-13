def avg_matrix_score(matrix, norm=False, cnorm=False):
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
