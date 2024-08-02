import pandas as pd
import sys
import logging

from settings import PYTHON_LOGGING_LEVEL, DIR_UKBB, DIR_UTILS, INTERM_FOLDER, IDP_RETINA_USED, IDP_VASCULAR_USED, Z_SCORE_APPLIED, COVARIATES, TYPE_BP, DATE_USED, MULTIPLE_TESTING, SQUARE_FIG, DOBLE_FIG, BP_USED, SOFT_MULTIPLE_TESTING, FILTER_OUTLIERS, DIR_OUTPUT, PLOT_SCATTER_FIGS, TYPE_OF_FIGS, SAVE_FIGURES

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)

sys.path.append(DIR_UTILS)

## IMPROTANT, if using MULTIPLE_TESTING = Subset the complete figure does not apply a proper multiple testing since I am interested in the sub set ()
# Since the image is Squared (IDPs, IDPs) + Rectangular (IDPs, Retina) -> Subset = IDPs * (IDPs/2 + Retina)

from utils.plotting.generate_plots import (
    filtered_corr_plots,
    pval_asterisk,
    doble_plot_corr_pval,
    plot_multiple_fig_shared,
    plot_corr_pval,
    filter_dataframe,
    pearsonr_pval,
)

from utils.plotting.plotting_settings import pheno_img_values
from utils.plotting.generate_plots import get_names_figs
from utils.preprocessing_pheno.compute_rsid import get_files_names
from utils.data_information.data_IDPs_info import *
from utils.data_information.data_homologous_info import *
from utils.data_information.data_retina_info import list_retina_homologous_red_new
from utils.plotting.generate_secundary_plots import scatter_plot_secondary_figs

SAVE_FIGURES = bool(SAVE_FIGURES)

if TYPE_OF_FIGS == "suppl_figs1":
    combine_IMT_heart_brain = False  # True # False -> plot the _corr_brain fig
    main_final_figs = False  # True

elif TYPE_OF_FIGS == "suppl_figs2":
    combine_IMT_heart_brain = True  # True # False -> plot the _corr_brain fig
    main_final_figs = False  # True

elif TYPE_OF_FIGS == "main_figs":
    combine_IMT_heart_brain = True  # True # False -> plot the _corr_brain fig
    main_final_figs = True  # True


save_dir = f"{DIR_UKBB}{INTERM_FOLDER}"


def plot_phenotypic_analysis():
    """
    This function performs phenotypic analysis by plotting correlation matrices and scatter plots
    for the given data. It reads data from CSV files, filters the data based on certain conditions,
    and generates various plots.

    Parameters:
        None

    Returns:
        None
    """
    name_file_values_count, name_file_zscored = get_files_names(IDP_RETINA_USED, 
                                                                IDP_VASCULAR_USED, 
                                                                Z_SCORE_APPLIED, 
                                                                combine_IMT_heart_brain, 
                                                                COVARIATES,
                                                                TYPE_BP,
                                                                DATE_USED,
                                                                FILTER_OUTLIERS)


    title_square, title_doble, title_doble_N = get_names_figs(DATE_USED, 
                                                SQUARE_FIG, 
                                                DOBLE_FIG, 
                                                Z_SCORE_APPLIED, 
                                                combine_IMT_heart_brain, 
                                                COVARIATES, 
                                                MULTIPLE_TESTING, 
                                                SAVE_FIGURES,
                                                TYPE_BP,
                                                FILTER_OUTLIERS,
                                                TYPE_OF_FIGS)

    df = pd.read_csv(f'{save_dir}{name_file_zscored}', index_col=False)
    df_non_null_values_count = pd.read_csv(f'{save_dir}{name_file_values_count}', index_col=False)


    df_non_null_values_count.index = df_non_null_values_count.columns

    if main_final_figs:
        filter_list_IDPs = list_retina_homologous_red_new + list(dict_final_IDPs.values())
        df = df[filter_list_IDPs]
        df_non_null_values_count = df_non_null_values_count.loc[filter_list_IDPs,filter_list_IDPs]

    if combine_IMT_heart_brain == False:

        dict_brain = dict_CBF.copy()
        dict_brain.update(dict_ATT)
        dict_brain.update(dict_deletions_brain)
        dict_brain.update(dict_brain_vessel)

        dict_heart = dict_heart_vessel.copy()
        dict_heart.update(dict_heart_vessel_2)
        dict_heart.update(dict_heart_functional)

        print(f"Plotting individual organ corr. Saving...")
        cte = 1.3
        filtered_corr_plots(
            df,
            dict_brain,
            figsize_1=(cte * 26, 26),
            cbar_1=0.05,
            SAVE_FIGURES=SAVE_FIGURES,
            name_used=f"Brain_cov_{COVARIATES}_outliers_{FILTER_OUTLIERS}",
            only_half=True,
        )
        filtered_corr_plots(
            df,
            dict_IMT,
            figsize_1=(cte * 10, 10),
            cbar_1=0.05,
            SAVE_FIGURES=SAVE_FIGURES,
            name_used=f"Carotid_cov_{COVARIATES}_outliers_{FILTER_OUTLIERS}",
            only_half=True,
        )
        filtered_corr_plots(
            df,
            dict_heart,
            figsize_1=(cte * 15, 15),
            cbar_1=0.05,
            SAVE_FIGURES=SAVE_FIGURES,
            name_used=f"Aorta_cov_{COVARIATES}_outliers_{FILTER_OUTLIERS}",
            only_half=True,
        )

    figsize_val_square, cbar_fraction_val_square, cte_square, figsize_val_both, width_ratios_val_both = pheno_img_values(IDP_VASCULAR_USED, combine_IMT_heart_brain, main_final_figs)
               
    list_values = [
        col
        for col in df.columns
        if col not in list_retina_homologous_red_new and col != "new_eid"
    ]

    #####
    if TYPE_OF_FIGS == "suppl_figs2":
        list_values = list(dict_homologous_averag.values())

    df = df[list_values + list_retina_homologous_red_new]

    df_corr_matrix = df.corr()
    plot_multiple_fig_shared(
        df_corr_matrix,
        df_non_null_values_count,
        figsize_val=figsize_val_square,
        cbar_fraction_val=cbar_fraction_val_square,
    )

    # ### P values plots
    # df.corr(method=spearmanr) and df.corr(method=pearsonr) needs to not have nans
    df_pval = df.corr(method=pearsonr_pval)

    ### Save csvs:
    df_corr_matrix.to_csv(f"{DIR_OUTPUT}{name_file_zscored[:-4]}_{TYPE_BP}_corr_matrix_{TYPE_OF_FIGS}.csv", index=False)
    df_pval.to_csv(f"{DIR_OUTPUT}{name_file_zscored[:-4]}_{TYPE_BP}_pval_matrix_{TYPE_OF_FIGS}.csv", index=False)

    df_val_asterisk = pval_asterisk(
        df,
        df_pval,
        MULTIPLE_TESTING,
        SOFT_MULTIPLE_TESTING,
        list_values,
        list_retina_homologous_red_new,
        log_pval=False,
        N_shape=False,
    )

    plot_corr_pval(df.corr(), df_val_asterisk, cte_square)

    red_df_corr_matrix = filter_dataframe(
        df_corr_matrix, list_values, list_retina_homologous_red_new
    )
    # red_df_non_null_values_count = filter_dataframe(df_non_null_values_count, list_retina_homologous_red_new, list_values)

    # red_df_corr_matrix = red_df_corr_matrix.T
    # red_df_non_null_values_count = red_df_non_null_values_count.T

    df_IDPs = df[list_values]
    df_IDPs_pval = df_IDPs.corr(method=pearsonr_pval)

    df_IDPs_pval_asterisk = pval_asterisk(
        df_IDPs,
        df_IDPs_pval,
        MULTIPLE_TESTING,
        SOFT_MULTIPLE_TESTING,
        list_values,
        list_retina_homologous_red_new,
        log_pval=False,
    )  # , N_shape=len(df_IDPs.columns)**2)

    # plot_corr_pval(df_IDPs.corr(),df_IDPs_pval_asterisk, val_conversion_factor=3, val_conversion_factor_extra=0.7, title_fig=title_square, only_half=True)

    red_df_val_asterisk = filter_dataframe(
        df_val_asterisk, list_values, list_retina_homologous_red_new
    )
    # red_df_val_asterisk = red_df_val_asterisk.T

    # plot_corr_pval(red_df_corr_matrix, red_df_val_asterisk, 5.5, val_conversion_factor_extra=0.5, title_fig=False, cmap_used='RdBu_r')

    doble_plot_corr_pval(
        df_corr1=df_IDPs.corr(),
        df_pval1=df_IDPs_pval_asterisk,
        df_corr2=red_df_corr_matrix,
        df_pval2=red_df_val_asterisk,
        figsize_val=figsize_val_both,
        width_ratios_val=width_ratios_val_both,
        title_fig=title_doble,
        only_half1=True,
        only_half2=False,
        cmap_used1="seismic",
        cmap_used2="RdBu_r",
        cbar_BP = TYPE_BP
    )

    if combine_IMT_heart_brain == True:
        red_df_non_null_values_count = df_non_null_values_count.loc[
            list_values, list_retina_homologous_red_new
        ]
        df_non_null_values_count_IDPs = df_non_null_values_count.loc[
            list_values, list_values
        ]
        # red_df_non_null_values_count = red_df_non_null_values_count.T
        doble_plot_corr_pval(
            df_corr1=df_non_null_values_count_IDPs,
            df_pval1=None,
            df_corr2=red_df_non_null_values_count,
            df_pval2=None,
            figsize_val=figsize_val_both,
            width_ratios_val=width_ratios_val_both,
            title_fig=title_doble_N,
            only_half1=True,
            only_half2=False,
            cmap_used1="Blues",
            cmap_used2="Greys"
        )  # Greys ()
        # cmax_used = df_non_null_values_count_IDPs.max().max())

    if PLOT_SCATTER_FIGS in ["True", "true"]:
        scatter_plot_secondary_figs(df)