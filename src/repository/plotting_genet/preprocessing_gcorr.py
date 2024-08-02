import sys
import logging
import numpy as np
from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS, DATE_USED, FILE_GENETIC_INFO, MULTIPLE_TESTING, SOFT_MULTIPLE_TESTING, DIR_OUTPUT

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)

sys.path.append(DIR_UTILS)

from utils.plotting.generate_plots import pval_asterisk
from utils.preprocessing_genet.previous_sumstats import (
    get_genetic_names,
    read_ldsr,
    get_gcorr_file_names,
    reorder_index_columns,
    rename_squared,
    squared_detele_col_index,
    convert_df_figures,
    detele_col_index,
)
from utils.data_information.data_genetics_info import list_names
from utils.data_information.data_retina_info import (
    list_retina_homologous_red,
    list_retina_homologous_red_new,
)

## TO DEFINE
complete_list = False  ## To have the squared figures or the subset
pseudo_homologous = True  # True # True

#dir_output_data_folder = f"{DIR_UKBB}{INTERM_FOLDER}"

## pheno_info_file is a csv to filter the relevant IDPs
pheno_info_file = f"{DIR_UTILS}utils/data_information/{FILE_GENETIC_INFO}"


def preprocesing_ldsr_output():
    """
    Preprocesses the output of LDSR analysis.

    This function performs various data manipulations and transformations on the output of LDSR analysis,
    including reordering index and columns, renaming dataframes, converting dataframes to CSV files, and
    performing statistical calculations.

    Returns:
        None
    """
    traits_all, traits_all_new = get_genetic_names(pheno_info_file)

    #logger.info(traits_all_new)

    # path This we can read from config
    df_diseases_all = get_gcorr_file_names(traits_all, list_retina_homologous_red)

    # datafields_irnt = [ dat + "_irnt.gwas.imputed_v3.both_sexes.tsv.sumstats.gz" for dat in traits_reduced]
    datafields_irnt = list(df_diseases_all["file"])
    datafields_pheno = [
        f"{dat}__munged.sumstats.gz" for dat in list_retina_homologous_red
    ]
    diseasess_tra_aux = list(df_diseases_all["pheno"])

    traits_col_index = list_retina_homologous_red + diseasess_tra_aux
    traits_names = datafields_pheno + datafields_irnt

    df_cov, df_corr, df_std2, df_pval, df_h2, df_h2_std = read_ldsr(traits_names, traits_col_index)

    ### Reorder index and columns
    # Apply the function to each DataFrame
    df_cov = reorder_index_columns(df_cov, list_names)
    df_corr = reorder_index_columns(df_corr, list_names)
    df_std2 = reorder_index_columns(df_std2, list_names)
    df_pval = reorder_index_columns(df_pval, list_names)
    df_h2 = reorder_index_columns(df_h2, list_names)
    df_h2_std = reorder_index_columns(df_h2_std, list_names)

    #save csv
    df_h2_std['index'] = df_h2_std.index
    #df_h2_std['columns'] = df_h2_std.columns
    df_h2_std['h2_std'] = np.diag(df_h2_std)
    df_h2_std[['index', 'h2_std']].to_csv(f"{DIR_OUTPUT}{DATE_USED}_h2_std.csv", index=False)
    df_corr.to_csv(f"{DIR_OUTPUT}{DATE_USED}_g_corr.csv", index=True)
    df_pval.to_csv(f"{DIR_OUTPUT}{DATE_USED}_g_pval.csv", index=True)



    # for squared plot:
    df_corr_all = rename_squared(
        df_corr,
        list_retina_homologous_red + traits_all,
        list_retina_homologous_red_new + traits_all_new,
    )
    df_std_all = rename_squared(
        df_std2,
        list_retina_homologous_red + traits_all,
        list_retina_homologous_red_new + traits_all_new,
    )
    df_pval_all = rename_squared(
        df_pval,
        list_retina_homologous_red + traits_all,
        list_retina_homologous_red_new + traits_all_new,
    )
    df_h2_all = rename_squared(
        df_h2,
        list_retina_homologous_red + traits_all,
        list_retina_homologous_red_new + traits_all_new,
    )

    df_corr_all_cont, df_corr_minus_std, df_std_all = convert_df_figures(
        df_corr_all, df_std_all
    )
    linear_log10p_copy3_cont = pval_asterisk(
        df_corr_all_cont,
        df_pval_all,
        MULTIPLE_TESTING,
        SOFT_MULTIPLE_TESTING,
        traits_all_new,
        list_retina_homologous_red_new,
    )

    df_corr_all_cont_T = df_corr_all_cont.T
    linear_log10p_copy3_cont_T = linear_log10p_copy3_cont.T
    df_h2_all_T = df_h2_all.T

    ##to csv
    df_corr_all_cont_T.to_csv(
        f"{DIR_OUTPUT}{DATE_USED}_complete_corr_gcorr.csv", index=True
    )
    linear_log10p_copy3_cont_T.to_csv(
        f"{DIR_OUTPUT}{DATE_USED}_complete_log10p_gcorr.csv", index=True
    )
    linear_log10p_copy3_cont_T.to_csv(
        f"{DIR_OUTPUT}{DATE_USED}_complete_log10p_gcorr.csv", index=True
    )
    df_h2_all_T.to_csv(
        f"{DIR_OUTPUT}{DATE_USED}_complete_h2_gcorr.csv", index=True
    )

    # for y IDPs x retina plot:
    df_corr_retina_IDPs = detele_col_index(df_corr, list_retina_homologous_red_new)
    df_std_retina_IDPs = detele_col_index(df_std2, list_retina_homologous_red_new)
    df_pval_retina_IDPs = detele_col_index(df_pval, list_retina_homologous_red_new)
    df_h2_retina_IDPs = detele_col_index(df_h2, list_retina_homologous_red_new)

    df_corr_retina_cont, df_corr_minus_std_retina_IDPs, df_std_retina_IDPs = (
        convert_df_figures(df_corr_retina_IDPs, df_std_retina_IDPs)
    )
    linear_log10p_copy3_retina_IDPs = pval_asterisk(
        df_corr_retina_cont,
        df_pval_retina_IDPs,
        MULTIPLE_TESTING,
        SOFT_MULTIPLE_TESTING,
        traits_all_new,
        list_retina_homologous_red_new,
    )

    # title_rectangular = f'{DATE_USED}_g_IDPs_retina_multest_{MULTIPLE_TESTING}'
    # f_.individual_star(df_corr_retina_cont.T, linear_log10p_copy3_retina_IDPs.T, 2.5, val_conversion_factor_extra=0.5, title_fig=title_rectangular, cmap_used='RdBu_r')

    # df_pval.to_csv('/SSD/home/sofia/gcorr_pval.csv')

    df_corr_IDPs = squared_detele_col_index(df_corr, traits_all_new)
    df_std_IDPs = squared_detele_col_index(df_std2, traits_all_new)
    df_pval_IDPs = squared_detele_col_index(df_pval, traits_all_new)
    df_h2_IDPs = squared_detele_col_index(df_h2, traits_all_new)

    df_corr_IDPs_cont, df_corr_minus_std_IDPs, df_std_IDPs = convert_df_figures(
        df_corr_IDPs, df_std_IDPs
    )
    linear_log10p_copy3_IDPs = pval_asterisk(
        df_corr_IDPs_cont,
        df_pval_IDPs,
        MULTIPLE_TESTING,
        SOFT_MULTIPLE_TESTING,
        traits_all_new,
        list_retina_homologous_red_new,
    )

    df_corr_IDPs_cont_T = df_corr_IDPs_cont.T
    linear_log10p_copy3_IDPs_T = linear_log10p_copy3_IDPs.T
    df_corr_retina_cont_T = df_corr_retina_cont.T
    linear_log10p_copy3_retina_IDPs_T = linear_log10p_copy3_retina_IDPs.T

    df_h2_IDPs_T = df_h2_IDPs.T
    df_h2_retina_IDPs_T = df_h2_IDPs.T
    # to csv
    df_corr_IDPs_cont_T.to_csv(
        f"{DIR_OUTPUT}{DATE_USED}_first_corr_gcorr.csv", index=True
    )
    linear_log10p_copy3_IDPs_T.to_csv(
        f"{DIR_OUTPUT}{DATE_USED}_first_log10p_gcorr.csv", index=True
    )
    df_corr_retina_cont_T.to_csv(
        f"{DIR_OUTPUT}{DATE_USED}_second_corr_gcorr.csv", index=True
    )
    linear_log10p_copy3_retina_IDPs_T.to_csv(
        f"{DIR_OUTPUT}{DATE_USED}_second_log10p_gcorr.csv", index=True
    )
    df_corr_retina_cont_T.to_csv(
        f"{DIR_OUTPUT}{DATE_USED}_second_corr_gcorr.csv", index=True
    )
    df_h2_IDPs_T.to_csv(
        f"{DIR_OUTPUT}{DATE_USED}_second_h2_gcorr.csv", index=True
    )
    df_h2_retina_IDPs_T.to_csv(
        f"{DIR_OUTPUT}{DATE_USED}_second_h2_gcorr.csv", index=True
    )
