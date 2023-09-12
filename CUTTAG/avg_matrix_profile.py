def avg_matrix_score(matrix, norm=False, cnorm=False):
    
    df = pd.read_csv(matrix, sep='\t', comment='@', header=None)
    df = df.drop(columns=list(range(6)))
    
    if norm:
        df[df>np.nanquantile(df, 0.99)] = np.nanquantile(df, 0.99)
    
    return df.mean(axis=0).values


