import pandas as pd
import sys
import logging

from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS, DATE_USED, INTERM_FOLDER, IDP_RETINA_USED, DIR_FILE_EID_MAP, DIR_FILE_RETINA_IDPS, OLD_EID, END_RAW_IDPS_COV_RETINA, DIR_UKBB

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)

sys.path.append(DIR_UTILS)

from utils.create_dataset.create_dataset import (
    show_file_and_datafields,
    check_phenotypes_values_count,
    merge_and_size,
    get_retina_pheno_names,
)
from utils.data_information.data_IDPs_info import (
    list_keys_brain,
    list_keys_heart,
    list_keys_carotid,
)
from utils.data_information.data_dataset_info import list_folders_files
from utils.data_information.data_homologous_info import dict_final_IDPs
from utils.data_information.data_covariants_info import (
    list_all_covar_keys,
    covariants_GPC,
)

#logger.info

list_idps_keys = list_keys_brain + list_keys_heart + list_keys_carotid
list_idps_keys_gion = [key + "-" for key in list_idps_keys]
list_keys = list_idps_keys + list_all_covar_keys

# file dir and name to save the final csv
dir_file_retina_idps_cov = (
    f"{DIR_UKBB}{INTERM_FOLDER}{DATE_USED}_{IDP_RETINA_USED}{END_RAW_IDPS_COV_RETINA}"
)

def create_dataset_csv():
    """
    Creates a dataset CSV file by selecting and merging data from multiple CSV files.

    Returns:
        None
    """
    # Go through list_folders_files checking where is the data of interest (list_keys)
    all_columns_found, columns_not_found = show_file_and_datafields(
        DIR_UKBB, list_folders_files, list_keys
    )
    #logger.info(f"Columns not found: {columns_not_found}")

    # Select the data of interest from each csv and merge them
    for index, file_uk in enumerate(list_folders_files):
        current_list = ["eid"] + all_columns_found[file_uk]
        ## GPCs:
        if file_uk == "dataset_674797/ukb674797.csv":
            current_list = current_list + list(covariants_GPC.keys())

        logger.info(f"Selecting data from {file_uk}")
        df_file = pd.read_csv(f"{str(DIR_UKBB)}{file_uk}", usecols=current_list)
        df_merge = df_file.copy() if index == 0 else df_merge.merge(df_file, on="eid") #inner by default
    unique_columns_count = len(df_merge.columns.unique())
    total_columns_count = len(df_merge.columns)
    #logger.info(f"len(df_merge.columns.unique())==len(df_merge.columns): {unique_columns_count == total_columns_count}")
    assert unique_columns_count == total_columns_count, f"Number of unique columns ({unique_columns_count}) does not match the total number of columns ({total_columns_count})"

    ### CHECK THE SAMPLE SIZE IS THE SAME AS IN THE UKBB WEBPAGE (this can be improved by reading from the table)
    logger.info("Sanity check:")
    check_phenotypes_values_count(df_merge)

    ### Retina
    list_retina, list_retina_new, default_eid, IDPs_retina_name = get_retina_pheno_names(
        IDP_RETINA_USED
    )

    df_eid = pd.read_csv(DIR_FILE_EID_MAP, header=None)
    df_eid.rename(columns={0: OLD_EID, 1: default_eid}, inplace=True)
    ### Load retina
    list_retina_df = [default_eid] + list_retina
    df_retina = pd.read_csv(DIR_FILE_RETINA_IDPS, usecols=list_retina_df)
    df_retina.rename(columns={default_eid: OLD_EID}, inplace=True)

    ### Merge: (idps and cov) with (retina)
    # NOTE IMPORTANT, here we follow the aproach of not filtering only the idps with retina data since we are interested in the cross organ
    df_merge_vascular = merge_and_size(df_eid, df_merge, default_eid, how_used="left")
    df_merge_eid_retina = merge_and_size(
        df_merge_vascular, df_retina, OLD_EID, how_used="left"
    )
    # rm old eid
    df_merge_eid_retina.drop(columns=[OLD_EID], inplace=True)

    name_idps_columns = [
        col
        for col in df_merge_eid_retina.columns
        if any(substring in col for substring in list_idps_keys_gion)
    ]
    logger.info(f"IDPs names: {name_idps_columns}")

    ## Check the number of subjects:
    logger.info("Check the num of subjects for the main ipds")
    main_idps = [key for key in dict_final_IDPs.keys() if key != "Min carotid IMT"]
    main_idps_2 = [f"{idp}-2.0" for idp in main_idps]
    for name in main_idps_2:
        count = df_merge_eid_retina[name].count()
        logger.info(name, count)

    logger.info(f"Saving values... {dir_file_retina_idps_cov}")
    df_merge_eid_retina.to_csv(dir_file_retina_idps_cov, index=False)