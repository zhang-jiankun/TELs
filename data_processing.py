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

##### insulation-related functions

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

##### Stripe-related functions 

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

##### Loop-related functions 

def overlap_multires_loop(bedpe_small, bedpe_large, slop=5000):
	_bedpe_small = bedpe_small.copy()
	_bedpe_large = bedpe_large.copy()
	_bedpe_small['name1'] = 'small_' + _bedpe_small.index.astype(str)
	_bedpe_small['name2'] = 'small_' + _bedpe_small.index.astype(str)
	_bedpe_large['name1'] = 'large_' + _bedpe_large.index.astype(str)
	_bedpe_large['name2'] = 'large_' + _bedpe_large.index.astype(str)
	bt_small = pbt.BedTool.from_dataframe(_bedpe_small)
	bt_large = pbt.BedTool.from_dataframe(_bedpe_large)
	overlap = bt_large.pairtopair(bt_small, type='both', slop=slop).to_dataframe(names=list(range(18)))
	assert (overlap[6] == overlap[7]).all() # filtered large 
	assert (overlap[14] == overlap[15]).all() # keeped small 
	filtered_large = np.unique(overlap[6].values)
	final_large = _bedpe_large.copy()
	final_large = final_large[~final_large['name1'].isin(filtered_large)]
	merged_bedpe = pd.concat([final_large, _bedpe_small], ignore_index=True)
	return merged_bedpe

def loop_classfication(chromloop, chromhmm='/dshare/xielab/analysisdata/ThreeDGene/zhangjk/GM12878-MicroC/proc_data/chromhmm/chromhmm_gm12878.bed', slop=5000): 
	"""
		chromhmm : ../proc_data/chromhmm/chromhmm_gm12878.bed
		chromloop: ../analyses/loop_calling/chromosight_loops/res_5000/chromosight_loops_res_5000.bedpe
	"""

	chromHmm = pd.read_csv(chromhmm, sep='\t', names=['chrom', 'start', 'end', 'state'])
	## chromHmm_ctcf = chromHmm[chromHmm['state']=='CTCF']
	chromHmm_tss  = chromHmm[chromHmm['state']=='TSS']
	chromHmm_enh  = chromHmm[(chromHmm['state']=='E') | (chromHmm['state']=='WE')]
	# bed_ctcf = pbt.BedTool.from_dataframe(chromHmm_ctcf)
	bed_tss  = pbt.BedTool.from_dataframe(chromHmm_tss)
	bed_enh  = pbt.BedTool.from_dataframe(chromHmm_enh)

	chromLoop = chromloop.copy()
	chromLoop['dot1'] = 'dot1-' + chromLoop.index.astype(str)
	chromLoop['dot2'] = 'dot2-' + chromLoop.index.astype(str)
	dot1_bed = pbt.BedTool.from_dataframe(chromLoop[['chr1', 'start1', 'end1', 'dot1']].sort_values(['chr1', 'start1'], ascending=True)).slop(b=slop, genome='hg38')
	dot2_bed = pbt.BedTool.from_dataframe(chromLoop[['chr2', 'start2', 'end2', 'dot2']].sort_values(['chr2', 'start2'], ascending=True)).slop(b=slop, genome='hg38')

	##### identify CTCF-CTCF loops
	# bin1_ctcf = dot1_bed.intersect(bed_ctcf, wa=True, u=True).to_dataframe()['name'].values 
	# bin2_ctcf = dot2_bed.intersect(bed_ctcf, wa=True, u=True).to_dataframe()['name'].values 
	# chromLoop_ctcf = chromLoop[(chromLoop.dot1.isin(bin1_ctcf)) & (chromLoop.dot2.isin(bin2_ctcf))][['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']]

	##### identify promoter-enhancer loops
	bin1_tss = dot1_bed.intersect(bed_tss, wa=True, u=True).to_dataframe()['name'].values
	bin2_enh = dot2_bed.intersect(bed_enh, wa=True, u=True).to_dataframe()['name'].values
	chromLoop_enhTss_1 = chromLoop[(chromLoop.dot1.isin(bin1_tss)) & (chromLoop.dot2.isin(bin2_enh))][['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']]
	bin1_enh = dot1_bed.intersect(bed_enh, wa=True, u=True).to_dataframe()['name'].values
	bin2_tss = dot2_bed.intersect(bed_tss, wa=True, u=True).to_dataframe()['name'].values
	chromLoop_enhTss_2 = chromLoop[(chromLoop.dot1.isin(bin1_enh)) & (chromLoop.dot2.isin(bin2_tss))][['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']]
	chromLoop_enhTss = pd.concat([chromLoop_enhTss_1, chromLoop_enhTss_2], ignore_index=True).sort_values(['chr1', 'start1', 'start2'], ascending=True).drop_duplicates(subset=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'], keep='first')
	chromLoop_enhTss_2u = pd.concat([chromLoop_enhTss, chromLoop_enhTss_1], ignore_index=True).sort_values(['chr1', 'start1', 'start2'], ascending=True).drop_duplicates(subset=['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'], keep=False)
	
	##### identify promoter-promoter loops
	chromLoop_tss = chromLoop[(chromLoop.dot1.isin(bin1_tss)) & (chromLoop.dot2.isin(bin2_tss))][['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']]

	return chromLoop_enhTss, chromLoop_enhTss_1, chromLoop_enhTss_2u, chromLoop_tss

def extract_ctcf_loop(chromloop):
	chromLoop = chromloop.copy()
	chromLoop['name'] = chromLoop.index
	dot1_bed = pbt.BedTool.from_dataframe(chromLoop[['chrom1', 'start1', 'end1', 'name']].sort_values(['chrom1', 'start1'], ascending=True))
	dot2_bed = pbt.BedTool.from_dataframe(chromLoop[['chrom2', 'start2', 'end2', 'name']].sort_values(['chrom2', 'start2'], ascending=True))

	ctcf_sites = pd.read_csv('~/database/JASPAR/results/GM12878/CTCF/ENCFF960ZGP.bed', sep='\t', header=None)[[0, 1, 2]]
	ctcf_sites.columns = ['chrom', 'start', 'end']
	pos_ctcf = pbt.BedTool.from_dataframe(ctcf_sites)

	bin1 = dot1_bed.intersect(pos_ctcf, wa=True, u=True).to_dataframe()['name'].values
	bin2 = dot2_bed.intersect(pos_ctcf, wa=True, u=True).to_dataframe()['name'].values
	return chromLoop[(chromLoop['name'].isin(bin1)) & (chromLoop['name'].isin(bin2))].drop(columns=['name'])
