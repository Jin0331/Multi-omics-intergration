"""
example : python src/Multi-omics-integration-analysis.py \
         -b /home/wmbio/WORK/gitworking/Multi-omics-intergration/ \
         -c BRCA

@author: Jinwoo Lee
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
from wmbio import * 
import argparse

if __name__ == "__main__": 

    parser = argparse.ArgumentParser(description='Subgroup Analysis!')
    parser.add_argument('-b', '--base', required=True, type=str, help='Root Path')
    parser.add_argument('-c', '--cancer', required=True, type=str, help='Types of cancer')
    parser.add_argument('-d', '--dea', default="deseq2", type=str, help='DESeq2(deseq2) or EdgeR(edger) or ALL(all)')
    parser.add_argument('-l', '--logfc', default=1, type=int, help='DESeq2(deseq2) or EdgeR(edger) or ALL(all)')
    parser.add_argument('-f', '--fdr', default=0.05, type=int, help='DESeq2(deseq2) or EdgeR(edger) or ALL(all)')

    static_args = parser.parse_args()

    # file path
    os.chdir(static_args.base)
    CANCER_TYPE = static_args.cancer
    METHOD = static_args.dea
    LOGFC = static_args.logfc
    FDR = static_args.fdr

    GROUP_PHTH = os.getcwd() + '/group/'
    PNG_PATH = os.getcwd() + '/png/'
    GROUP_VALIDATION_PATH = os.getcwd() + '/group_validation/'
    DEG_PATH = os.getcwd() + "/best_deg/"
    RDATA_PATH = os.getcwd() + "/RAW_DATA/GDC_PREPARE/"


    # Load Validation score
    col=['FILENAME','Log Rank Test','Silhouette','RNA_ANOVA_F1','RNA_RF_F1',
        'miRNA_ANOVA_F1','miRNA_RF_F1','Methylation_ANOVA_F1','Methylation_RF_F1']

    group_score = pd.read_csv(GROUP_VALIDATION_PATH + CANCER_TYPE + '_validation.csv', usecols=col)

    # Condition for Filtering
    # 
    filter_cond = (group_score['Silhouette'] >= 0.2) & (group_score['Log Rank Test'] < 0.05) & \
                ((group_score['RNA_ANOVA_F1'] > 80) | (group_score['RNA_RF_F1'] > 80)) & \
                ((group_score['miRNA_ANOVA_F1'] > 80) | (group_score['miRNA_RF_F1'] > 80)) & \
                ((group_score['Methylation_ANOVA_F1'] > 85) | (group_score['Methylation_RF_F1'] > 85))
    group_score = group_score[filter_cond].sort_values(["Silhouette"], ascending = (False))

    bestSubgroup = group_score.FILENAME.to_list()

    # 5 fold로 제한.
    if len(bestSubgroup) >= 100:
        bestSubgroup = bestSubgroup[:100]
    print("SubGroup count : ", len(bestSubgroup))

    # DEA
    dea_result = list()
    for best_group in bestSubgroup:
        
        DEG_CHECK = "_".join([CANCER_TYPE, METHOD.upper(), best_group]) + ".txt"
        SAMPLE_GROUP = GROUP_PHTH + CANCER_TYPE + "/" + CANCER_TYPE + "_GROUP_" + best_group + ".txt"
        
        if os.path.isfile(DEG_PATH + CANCER_TYPE + "/" + DEG_CHECK):
            deg_list = pd.read_csv(DEG_PATH + CANCER_TYPE + "/" + DEG_CHECK, sep = "\t")
        else :
            # DEG Extraction
            deg_list = deg_extract(log_fc=LOGFC, fdr=FDR,
                          cancer_type=CANCER_TYPE, 
                          sample_group=SAMPLE_GROUP, deg_path=DEG_PATH, 
                          file_name=best_group,
                          rdata_path=RDATA_PATH,
                          method=METHOD,
                          batch_removal=True)
        
        dea_result.append(deg_list)

    # Filter DEA
    # combine result
    if METHOD == 'all':
        dea_combine = list(map(deseq2_edger_combine, dea_result))
        dea_combine = [col_rename(dea_combine[index], index, bestSubgroup) for index in range(len(dea_combine))]
        dea_combine = reduce(lambda left, right : pd.merge(left, right, left_on='gene', right_on='gene', how = 'outer'), dea_combine)
    elif METHOD == 'deseq2' :
        dea_combine = list(map(lambda d : d[["row", "log2FoldChange"]], dea_result))
        dea_combine = [col_rename(dea_combine[index], index, bestSubgroup) for index in range(len(dea_combine))]
        dea_combine = reduce(lambda left, right : pd.merge(left, right, left_on='gene', right_on='gene', how = 'inner'), dea_combine)

    # median & mean
    dea_combine["SubGroup-log2FC_median"] = dea_combine.iloc[:, 1:].median(axis=1)
    dea_combine["SubGroup-log2FC_mean"] = dea_combine.iloc[:, 1:].mean(axis=1)


    # NT vs TP DEA
    # log fc, FDR value 수정
    nt_tp_deseq2 = deg_extract_normal(log_fc=0, pvalue=0.1, cancer_type=CANCER_TYPE, 
                                  rdata_path=RDATA_PATH, deg_path=DEG_PATH, batch_removal=True)

    nt_tp_deseq2_col = nt_tp_deseq2[['row', 'log2FoldChange', 'pvalue']]
    nt_tp_deseq2_col.columns = ['gene', 'NT-TP_log2FoldChange', 'pvalue']

    result_combine = pd.merge(left=dea_combine, right=nt_tp_deseq2_col, left_on='gene', right_on='gene', how = 'left')   

    # textmining
    sql = 'SELECT * FROM ' + CANCER_TYPE
    tm_df = query_tm_db(sql)
    result_combine_tm = pd.merge(left=result_combine, right=tm_df, left_on="gene", right_on="gene", how='left')

    # Result write
    Path(os.getcwd() + "/RESULT").mkdir(parents=True, exist_ok=True)
    Path(os.getcwd() + "/RESULT/" + CANCER_TYPE).mkdir(parents=True, exist_ok=True)
    time_stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

    result_combine_tm.to_csv(os.getcwd() + "/RESULT/" + CANCER_TYPE + "/" + CANCER_TYPE + '-' + time_stamp +'.csv', index = False)
