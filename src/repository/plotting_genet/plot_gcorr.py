import sys
import logging
import numpy as np
from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS, DATE_USED, MULTIPLE_TESTING, SOFT_MULTIPLE_TESTING, DIR_OUTPUT

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)

sys.path.append(DIR_UTILS)

from utils.preprocessing_genet.previous_sumstats import (
    individual_star,
    read_and_process_csv,
    read_and_process_log10p_csv,
    doble_plot_corr_pval,
)
from utils.plotting.plotting_settings import FIG_PHENO_GENO_SIZE, FIG_PHENO_GENO_RATIOS


def plot_ldsr_genetic_corr():
    """
    Plots the genetic correlation between different datasets.

    Reads and processes CSV files containing correlation and p-value data for different datasets.
    Generates a double plot showing the correlation and p-values for IDPs and retina datasets.
    Adds individual stars to indicate statistical significance.

    Args:
        None

    Returns:
        None
    """
    
    # Read and process the CSV files
    df_corr_all = read_and_process_csv(
        f"{DIR_OUTPUT}{DATE_USED}_complete_corr_gcorr.csv"
    )
    df_log10p_all = read_and_process_log10p_csv(
        f"{DIR_OUTPUT}{DATE_USED}_complete_log10p_gcorr.csv"
    )

    df_corr_IDPs = read_and_process_csv(
        f"{DIR_OUTPUT}{DATE_USED}_first_corr_gcorr.csv"
    )
    df_log10p_IDPs = read_and_process_log10p_csv(
        f"{DIR_OUTPUT}{DATE_USED}_first_log10p_gcorr.csv"
    )

    df_corr_retina = read_and_process_csv(
        f"{DIR_OUTPUT}{DATE_USED}_second_corr_gcorr.csv"
    )
    df_log10p_retina = read_and_process_log10p_csv(
        f"{DIR_OUTPUT}{DATE_USED}_second_log10p_gcorr.csv"
    )
    heritability_diag=True
    if heritability_diag:
        #df_h2_all= read_and_process_csv(f'{DIR_OUTPUT}{DATE_USED}_complete_h2_gcorr.csv')
        df_h2_IDPs = read_and_process_csv(f'{DIR_OUTPUT}{DATE_USED}_second_h2_gcorr.csv')
        #df_h2_retina= read_and_process_csv(f'{DIR_OUTPUT}{DATE_USED}_second_h2_gcorr.csv')

        #Replace the diag of corr idps by h2
        diagonal_h2_IDPs = np.diag(df_h2_IDPs)
        np.fill_diagonal(df_corr_IDPs.values, diagonal_h2_IDPs)
        
        #Replace the diag of log10p idps by ''
        np.fill_diagonal(df_log10p_IDPs.values, '')

    title_doble = f"{DATE_USED}_g_IDPs_IDPs_and_IDPs_retina_multest_{MULTIPLE_TESTING}.jpg"
    doble_plot_corr_pval(
        df_corr1=df_corr_IDPs,
        df_pval1=df_log10p_IDPs,
        df_corr2=df_corr_retina,
        df_pval2=df_log10p_retina,
        figsize_val=FIG_PHENO_GENO_SIZE,  # (12, 4), #(14, 6)
        width_ratios_val=FIG_PHENO_GENO_RATIOS,  # [1,1], #[1, 0.7],
        title_fig=title_doble,
        only_half1=True,
        only_half2=False,
        cmap_used1="seismic",
        cmap_used2="RdBu_r",
    )

    individual_star(df_corr_all, df_log10p_all, 1.5)
