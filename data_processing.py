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

def filter_spenn(tmp, tech, res=5000):
    
    """
        chr, pos1, pos2: anchor;
        chr2, pos3, pos4: span;
        notebooks/Suppl_Figure_Cool_analysis.ipynb
    """
    
    tmp.pos1 = tmp.pos1 - 1
    tmp.pos3 = tmp.pos3 - 1
    tmp['score'] = '.'
    
    tmp_left = tmp[tmp.pos1 == tmp.pos3]
    tmp_right = tmp[tmp.pos2 == tmp.pos4]
    
    tmp_left['strand'] = '+'
    tmp_left = tmp_left.sort_values(['chr', 'pos1'], ascending=True)
    tmp_left = tmp_left.reset_index().drop(columns='index')
    tmp_left['name'] = f'{tech}_Stripenn_ha_' + tmp_left.index.astype(str)
    
    tmp_right['strand'] = '-'
    tmp_right = tmp_right.sort_values(['chr2', 'pos3'], ascending=True)
    tmp_right = tmp_right.reset_index().drop(columns='index')
    tmp_right['name'] = f'{tech}_Stripenn_va_' + tmp_right.index.astype(str)
    
    tmp_anchors = pd.concat([
        tmp_left[['chr', 'pos1', 'pos2', 'name', 'score', 'strand']],
        tmp_right[['chr', 'pos1', 'pos2', 'name', 'score', 'strand']],
    ], ignore_index=True).sort_values(['chr', 'pos1'], ascending=True).rename(columns={'chr': 'chr', 'pos1': 'start', 'pos2': 'end'})
    
    tmp_span = pd.concat([
        tmp_left[['chr2', 'pos3', 'pos4', 'name', 'score', 'strand']],
        tmp_right[['chr2', 'pos3', 'pos4', 'name', 'score', 'strand']]
    ], ignore_index=True).sort_values(['chr2', 'pos3'], ascending=True).rename(columns={'chr2': 'chr', 'pos3': 'start', 'pos4': 'end'})
    
    return tmp_anchors, tmp_span

def filter_specall(tmp, tech, res=5000):
    
    tmp_left = tmp[tmp.end1 - tmp.start1 == res].copy(deep=True)
    tmp_left = tmp_left.sort_values(['chr1', 'start1', 'end1', 'end2'], ascending=False)
    tmp_left = tmp_left.drop_duplicates(subset=['chr1', 'start1', 'end1'], keep='first')
    tmp_left = tmp_left[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']]
    tmp_left = tmp_left[tmp_left.end2 - tmp_left.start1 < 2000000]
    tmp_left = tmp_left.reset_index().drop(columns='index')
    
    tmp_left['score'] = '.'
    tmp_left['name'] = f'{tech}_Stripecaller_ha_' + tmp_left.index.astype(str)
    tmp_left['strand'] = '+'

    tmp_right = tmp[tmp.end2 - tmp.start2 == res].copy(deep=True)
    tmp_right = tmp_right.sort_values(['chr2', 'start2', 'end2', 'start1'], ascending=True)
    tmp_right = tmp_right.drop_duplicates(subset=['chr2', 'start2', 'end2'], keep='first')
    tmp_right = tmp_right[['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']]
    tmp_right = tmp_right[tmp_right.end2 - tmp_right.start1 < 2000000]
    tmp_right = tmp_right.reset_index().drop(columns='index')
    
    tmp_right['score'] = '.'
    tmp_right['name'] = f'{tech}_Stripecaller_va_' + tmp_right.index.astype(str)
    tmp_right['strand'] = '-'

    tmp_anchors = pd.concat([
        tmp_left[['chr1', 'start1', 'end1', 'name', 'score', 'strand']].rename(columns={'chr1': 'chr', 'start1': 'start', 'end1': 'end'}),
        tmp_right[['chr2', 'start2', 'end2', 'name', 'score', 'strand']].rename(columns={'chr2': 'chr', 'start2': 'start', 'end2': 'end'}),
    ], ignore_index=True).sort_values(['chr', 'start'], ascending=True)

    tmp_span = pd.concat([
        tmp_left[['chr1', 'start1', 'end2', 'name', 'score', 'strand']],
        tmp_right[['chr1', 'start1', 'end2', 'name', 'score', 'strand']]
    ], ignore_index=True).rename(columns={'chr1': 'chr', 'start1': 'start', 'end2': 'end'}).sort_values(['chr', 'start'], ascending=True)
    
    return tmp_anchors, tmp_span

