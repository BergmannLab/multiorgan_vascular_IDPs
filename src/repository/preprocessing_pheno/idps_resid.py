import itertools
import pandas as pd
#import numpy as np
import os
import sys
import logging

from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS, DIR_UKBB, INTERM_FOLDER, DATE_USED, IDP_RETINA_USED, IDP_VASCULAR_USED, Z_SCORE_APPLIED, BP_USED, TYPE_BP, FILTER_OUTLIERS, PLOT_HISTOGRAMS, HIST_FILE, NUM_STD
from utils.data_information.data_homologous_info import dict_final_IDPs

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)


if BP_USED in {"True", True}:
    FILE_RAW_COVARIANTS = os.getenv("FILE_RAW_COVARIANTS_BP")
    l_covariates = [False, TYPE_BP]
else:
    FILE_RAW_COVARIANTS = os.getenv("FILE_RAW_COVARIANTS")
    l_covariates = [False, "all"]  # , 'age_age2_sex'] # 'age_age2_sex', age_age2_sex_PCs, all

# os.chdir(DIR_UTILS) # for jupyternotebook
sys.path.append(DIR_UTILS)  # for script
from utils.preprocessing_pheno.preprocessing_pheno import (
    get_phenotypes_names,
    create_counts_df,
    replace_outliers_with_nan,
    replace_outliers_IQR_with_nan,
    get_baselines
)
from utils.preprocessing_pheno.compute_rsid import (
    get_files_names,
    get_cov_names,
    df_IDPs_cov_substracted,
)
from utils.plotting.generate_plots import plot_cov_histograms

NUM_STD = int(NUM_STD)

l_combine = [False, True]  # l_combine = [True]

save_dir = f"{DIR_UKBB}{INTERM_FOLDER}"

saved_name_file = f"{DATE_USED}_{IDP_VASCULAR_USED}_{IDP_RETINA_USED}.csv"
dir_saved_name_file = f"{save_dir}{saved_name_file}"

def calculate_resid_from_idps_csv():  # sourcery skip: low-code-quality
    """
    Calculate residuals from IDPs CSV file.

    This function reads a CSV file containing IDPs (Independent Variables) data, performs data preprocessing steps,
    and calculates residuals. The residuals are then saved to a new CSV file.

    Returns:
        None
    """
    for covar, combine_IMT_heart_brain in itertools.product(
        l_covariates, l_combine
    ):
        name_file_values_count, name_file_zscored = get_files_names(
            IDP_RETINA_USED,
            IDP_VASCULAR_USED,
            Z_SCORE_APPLIED,
            combine_IMT_heart_brain,
            covar,
            False,
            DATE_USED,
            FILTER_OUTLIERS
        )
        # if name_file_zscored=='2024_02_19_homologous_all_all_zscored_True_combine_True_cov_used_all.csv':
        #     logger.info(name_file_zscored)
        ### FILTER
        (
            list_retina,
            list_retina_new,
            eid,
            IDPs_retina_name,
            list_keys,
            list_values,
        ) = get_phenotypes_names(
            IDP_RETINA_USED, IDP_VASCULAR_USED, combine_IMT_heart_brain
        )
        # logger.info(list_values)

        df = pd.read_csv(dir_saved_name_file)

        column_counts = []
        for col in df.columns:
            count = df[col].count()
            column_counts.append((col, count))
            logger.info(col, count)
            if count < 3000:
                logger.error(
                    f"Error: The column {col} has fewer than 3000 non-null values ({count})"
                )
        ## option to filter the subjects with more than 10std from the main:
        if FILTER_OUTLIERS not in {"False", False}:
            if FILTER_OUTLIERS == "std_method":
                #10 std was already removed for the retina
                df = replace_outliers_with_nan(df, list_values, NUM_STD, save_dir)
                #df = replace_outliers_with_nan(df, list_values+list_retina_new, NUM_STD, save_dir)
            elif FILTER_OUTLIERS == "iqr_method":
                df = replace_outliers_IQR_with_nan(df, list_values+list_retina_new, save_dir)
            else:
                logger.error("Error: Invalid value for FILTER_OUTLIERS. Valid values are 'std_method' and 'iqr_method'.")

        #### Plot histograms:
        if PLOT_HISTOGRAMS in ["True", True]:
            name_hist = f"{DATE_USED}_IDPs_filter_outliers_{FILTER_OUTLIERS}_{HIST_FILE}"
            plot_cov_histograms(df[list_values+list_retina_new], name_hist)

        if (covar == False) or (covar == "False"):
            #### Normalize df
            if Z_SCORE_APPLIED:
                df = df.apply(
                    lambda x: (
                        (x - x.mean()) / x.std()
                        if x.notna().any() and x.name != eid
                        else x
                    ),
                    axis=0,
                )
                # df = df.apply(lambda x: (x - x.mean()) / x.std() if x.notna().any() else x, axis=0)
        else:

            list_cov_columns = get_cov_names(covar)

            # If if ((covar!=False) and (covar!='all')): regress out that cov from retina
            # if covar not in [False, "all"]:
            #     df = df_IDPs_cov_substracted(df, covar, list_retina_new)
            #     list_cov_columns = list_cov_columns + [covar]

            # filter
            df = df[[eid] + list_retina_new + list_values + list_cov_columns]

            if combine_IMT_heart_brain in [True, 'True']:
                df_baselines = get_baselines(df, list(dict_final_IDPs.values()))
                print(df_baselines)
                
            #### Normalize df
            if Z_SCORE_APPLIED:
                df = df.apply(
                    lambda x: (
                        (x - x.mean()) / x.std()
                        if x.notna().any() and x.name != eid
                        else x
                    ),
                    axis=0,
                )
            
            # Only from IDPs, retina it is alreay withouth them
            df = df_IDPs_cov_substracted(df, list_cov_columns, list_values)


        # df.count().sort_values(ascending=False)  # create_counts_df(df)

        #### Filter by eid, retina and idps (since we do not want the cov data anymore):
        df = df[[eid] + list_retina_new + list_values]

        # Compute the non null values intersccions
        df_non_null_values_count = create_counts_df(df[list_retina_new + list_values])

        df.to_csv(f"{save_dir}{name_file_zscored}", index=False)
        df_non_null_values_count.to_csv(
            f"{save_dir}{name_file_values_count}", index=False
        )
        logger.info(
            f"Files saved! {save_dir}{name_file_zscored} \n combine_IMT_heart_brain={combine_IMT_heart_brain}"
        )
