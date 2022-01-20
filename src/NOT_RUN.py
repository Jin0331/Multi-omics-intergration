def load_preprocess_tcga_dataset(pkl_path, raw_path, group_path, cancer_type, norm):
    
    if os.path.isfile(pkl_path + "/" + cancer_type + "_omics.pkl"):
        # sep
        omics = pd.read_pickle(pkl_path + "/" + cancer_type + "_omics.pkl")
        rna = pd.read_pickle(pkl_path + "/" + cancer_type + "_rna.pkl")
        mirna = pd.read_pickle(pkl_path + "/" + cancer_type + "_mirna.pkl")
        mt_join_gene_filter = pd.read_pickle(pkl_path + "/" + cancer_type + "_mt.pkl")
        
    else :
        # RNA gene expression
        col = pd.read_csv(raw_path + "tcga_RSEM_Hugo_norm_count.gz",
                     sep = "\t", index_col=0, nrows=0).columns.to_list()
        use_col = ['sample'] + cancer_select(cols=col, cancer_type=cancer_type)
        df_chunk = pd.read_csv(raw_path + "tcga_RSEM_Hugo_norm_count.gz",
                     sep = "\t", index_col=0, iterator=True, chunksize=50000, usecols=use_col)
        rna = pd.concat([chunk for chunk in df_chunk])
        rna = rna[rna.index.isin(non_zero_column(rna))].T
        
        rna.to_pickle(pkl_path + "/" + cancer_type + "_rna.pkl")

        # miRNA expression
        col = pd.read_csv(raw_path + "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz",
                     sep = "\t", index_col=0, nrows=0).columns.to_list()
        use_col = ['sample'] + cancer_select(cols=col, cancer_type=cancer_type)

        df_chunk = pd.read_csv(raw_path + "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz",
                         sep = "\t", index_col=0, iterator=True, chunksize=50000, usecols=use_col)
        mirna = pd.concat([chunk for chunk in df_chunk])
        mirna = mirna[mirna.index.isin(non_zero_column(mirna))].T
        
        mirna.to_pickle(pkl_path + "/" + cancer_type + "_mirna.pkl")

        # methylation
        col = pd.read_csv(raw_path + "jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz",
                     sep = "\t", index_col=0, nrows=0).columns.to_list()
        use_col = ['sample'] + cancer_select(cols=col, cancer_type=cancer_type)

        df_chunk = pd.read_csv(raw_path + "jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz",
                         sep = "\t", index_col=0, iterator=True, chunksize=50000, usecols=use_col)
        mt = pd.concat([chunk for chunk in df_chunk])

        mt_map = pd.read_csv(raw_path + "probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy", sep="\t")

        mt_join = pd.merge(mt, mt_map, how = "left", left_on = "sample", right_on = "#id")\
                 .drop(['chrom', 'chromStart', 'chromEnd', 'strand', '#id'], axis=1)
        mt_join = mt_join[mt_join.gene != "."]
        mt_join.dropna(subset = ["gene"], inplace=True)

        # gene mean 
        mt_join_gene_filter = mt_join.groupby(['gene']).mean()
        mt_join_gene_filter = mt_join_gene_filter[mt_join_gene_filter.index.isin(non_zero_column(mt_join_gene_filter))].T
        
        mt_join_gene_filter.to_pickle(pkl_path + "/" + cancer_type + "_mt.pkl")
        
        # set same column for merge
        rna['sample'] = rna.index
        mirna['sample'] = mirna.index
        mt_join_gene_filter['sample'] = mt_join_gene_filter.index
        
    # set column for unique
    rna.columns = list(map(lambda col : col + "_RNA", rna.columns.to_list()))
    mirna.columns = list(map(lambda col : col + "_miRNA", mirna.columns.to_list()))
    mt_join_gene_filter.columns = list(map(lambda col : col + "_Methylation", mt_join_gene_filter.columns.to_list()))
        
    # set index
    rna_index = rna.index.to_list()
    mirna_index = mirna.index.to_list()
    mt_join_gene_filter_index = mt_join_gene_filter.index.to_list()
    
    # normalization
    if norm:
        scalerX = StandardScaler()
        rna_scale = scalerX.fit_transform(rna)
        mirna_scale = scalerX.fit_transform(mirna)
        mt_join_gene_filter_scale = scalerX.fit_transform(mt_join_gene_filter)

        # missing impute
        imputer = KNNImputer(n_neighbors=10)
        rna_impute = imputer.fit_transform(rna_scale)
        mirna_impute = imputer.fit_transform(mirna_scale)
        mt_join_gene_filter_impute = imputer.fit_transform(mt_join_gene_filter_scale)

        # Pandas
        rna = pd.DataFrame(rna_impute, columns=rna.columns)
        rna.index = rna_index

        mirna = pd.DataFrame(mirna_impute, columns=mirna.columns)
        mirna.index = mirna_index

        mt = pd.DataFrame(mt_join_gene_filter_impute, columns=mt_join_gene_filter.columns)
        mt.index = mt_join_gene_filter_index
    else :
        # missing impute
        imputer = KNNImputer(n_neighbors=10)
        rna_impute = imputer.fit_transform(rna)
        mirna_impute = imputer.fit_transform(mirna)
        mt_join_gene_filter_impute = imputer.fit_transform(mt_join_gene_filtere)

        # Pandas
        rna = pd.DataFrame(rna_impute, columns=rna.columns)
        rna.index = rna_index

        mirna = pd.DataFrame(mirna_impute, columns=mirna.columns)
        mirna.index = mirna_index

        mt = pd.DataFrame(mt_join_gene_filter_impute, columns=mt_join_gene_filter.columns)
        mt.index = mt_join_gene_filter_index
    
    # omics join group
    group_df = pd.read_csv(group_path, sep = "\t", index_col=0)
    
    omics = [rna, mirna, mt]
    omics = list(map(lambda df : pd.merge(left=group_df, right=df, how="inner", 
                                          left_index=True, right_index=True), omics))
    
    # list to dict
    zipbObj = zip(data_type, omics)
    omics = dict(zipbObj)
    
    return omics




intersect_feature = list()

intersect_feature.append(list(set(feature_result["rna"][0][0].iloc[:,1].to_list()) & 
                              set(feature_result["rna"][0][0].iloc[:,1].to_list())))
intersect_feature.append(list(set(feature_result["mirna"][0][0].iloc[:,1].to_list()) & 
                              set(feature_result["mirna"][0][0].iloc[:,1].to_list())))
intersect_feature.append(list(set(feature_result["mt"][0][0].iloc[:,1].to_list()) & 
                              set(feature_result["mt"][0][0].iloc[:,1].to_list())))

intersect_feature_list = ["group","OS","OS.time"] + list(chain(*intersect_feature))


