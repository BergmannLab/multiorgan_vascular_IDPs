import pandas as pd
import sys
import logging

from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS, DIR_UKBB, INTERM_FOLDER, DATE_USED, END_RAW_IDPS_COV_RETINA, IDP_RETINA_USED, IDP_VASCULAR_USED   

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)

sys.path.append(DIR_UTILS)

from utils.preprocessing_pheno.preprocessing_pheno import (
    get_retina_names,
    get_IDPs_names,
    compare_instances_2_3,
    check_instances,
    combine_IMT_cols,
    combine_CBF_ATT_cols,
    covariants_preprocessing,
    check_and_filter_columns,
    rename_columns,
    check_final_columns,
    adding_extra_cov,
    check_column_counts,
    check_column_thresholds,
    #combine_instances_23,
    #check_column_counts,
    check_keys,
    check_vals,
)
from utils.create_dataset.create_dataset import merge_and_size
from utils.data_information.data_covariants_info import (
    list_all_covar_keys,
    list_all_covar_values,
    covariants_GPC,
)
from utils.data_information.data_homologous_info import dict_final_IDPs


save_dir = f"{DIR_UKBB}{INTERM_FOLDER}"
dir_file_retina_idps_cov = (
    f"{save_dir}{DATE_USED}_{IDP_RETINA_USED}{END_RAW_IDPS_COV_RETINA}"
)
# Read file with the eid correspondance, file with vascular IDPs, and retina IDPs:
SAVED_NAME_FILE = f"{DATE_USED}_{IDP_VASCULAR_USED}_{IDP_RETINA_USED}.csv"
DIR_SAVED_NAME_FILE = f"{save_dir}{SAVED_NAME_FILE}"
DIR_SAVED_NAME_STAS = f"{save_dir}{DATE_USED}_statistics_raw.csv"


def preproces_idps_cov():
    list_retina, list_retina_new, eid, IDPs_retina_name = get_retina_names(
        IDP_RETINA_USED
    )

    list_idp_keys, list_idp_values = get_IDPs_names(IDP_VASCULAR_USED, combine=False)
    list_keys = list_idp_keys + list_all_covar_keys
    list_values = list_idp_values + list_all_covar_values

    #### Load data: idps, retina and cov
    df_init = pd.read_csv(dir_file_retina_idps_cov)

    # covars preprocessing
    df = covariants_preprocessing(df_init)

    # idps preprocessing
    ## Get the complete idp names
    name_idps_columns = [
        col
        for col in df.columns
        if any(substring in col for substring in list_idp_keys)
    ]
    print(f"Num of cols: {len(name_idps_columns)}")#, with names: {name_idps_columns}")

    df_vascular = df[[eid] + name_idps_columns]

    ### check that df_vascular only has 2.0 and 3.0:
    check_instances(df_vascular)
    # getting basic info about idps instances 2.0 and 3.0
    compare_instances_2_3(df_vascular, DIR_SAVED_NAME_STAS)

    ### NOT USED (COMBINE INSTANCES 2 AND 3) Merge the phentypes with instances 2.0 and 3.0 into one. 2.0, 3.0 and the average if both
    # df_vascular, param_names = combine_instances_23(df_vascular) #.combine_instances_23(df_vascular)
    # check_column_counts(df_vascular, '22679-2.0', ['22679'], min_count=48524)
    # check_column_counts(df_vascular, '22420-2.0', ['22420'])
    # check_column_counts(df_vascular, '24383', ['24383-2.0', '24383-3.0'])

    # Select only the columns that end with '-2.0' or contain 'eid'
    df_vascular_2 = df_vascular.filter(regex="-2.0$|" + eid)
    print(
        f"By filtering instance 2, the number of idps went from {len(df_vascular.columns)} to {len(df_vascular_2.columns)}"
    )

    # Rename columns by removing '-2.0'
    new_columns = [col.replace("-2.0", "") for col in df_vascular_2.columns]
    df_vascular_2.columns = new_columns

    # Merge idps, cov and retina back
    df_filtered = merge_and_size(df, df_vascular_2, eid, how_used="left")

    ## Check the number of subjects:
    check_keys(df_filtered, dict_final_IDPs)

    # Rename idps, cov and retina cols, and check the cols
    try:
        df_merge_eid_retina_red = check_and_filter_columns(
            df_filtered, list_retina, list_keys, list_values, eid
        )
        check_keys(df_merge_eid_retina_red, dict_final_IDPs)
        if df_merge_eid_retina_red is not None:
            df_merge_eid_retina_red = rename_columns(
                df_merge_eid_retina_red,
                list_keys,
                list_values,
                list_retina,
                list_retina_new,
            )
            check_vals(df_merge_eid_retina_red, dict_final_IDPs)
            check_final_columns(
                df_merge_eid_retina_red, list_retina_new, list_values, eid
            )
    except Exception as e:
        print("An error occurred:", e)

    # filter by retina
    df_merge_eid_retina_red = df_merge_eid_retina_red[
        [eid] + list_retina_new + list_values
    ]

    # Define the column thresholds
    column_thresholds = {
        "Min carotid IMT at 240 degrees": 13189,
        "Asc aorta min area": 45533
    }
    # Check the columns against their thresholds
    try:
        check_column_thresholds(df_merge_eid_retina_red, column_thresholds)
    except ValueError as e:
        print(e)

    ## Average IMT (across the multiple angles)
    df_merge_eid_retina_red_avg = combine_IMT_cols(df_merge_eid_retina_red)

    ## Average CBF and ATT (across the multiple brain parts)
    if IDP_VASCULAR_USED == "all":
        df_merge_eid_retina_red_avg = combine_CBF_ATT_cols(df_merge_eid_retina_red_avg)

    # Adding age**2, hypertense, etc
    df_merge_eid_retina_red_avg = adding_extra_cov(df_merge_eid_retina_red_avg)

    # Check the non-null counts and collect the data
    check_column_counts(df_merge_eid_retina_red_avg, min_count=3000)

    # Create a DataFrame from the collected data
    # df_column_counts = pd.DataFrame(column_counts, columns=['Column', 'Non-Null Count'])

    # Save the DataFrame to a CSV file
    # df_column_counts.to_csv('column_non_null_counts.csv', index=False)
    # print("CSV file with non-null counts has been saved successfully.")

    ## format
    # Ensure 'eid' is preserved and other columns are converted to float
    eid_column = df_merge_eid_retina_red_avg["eid"]
    other_columns = df_merge_eid_retina_red_avg.drop(columns=["eid"]).astype(float)
    # Combine 'eid' column back with the float-converted data
    df_merge_eid_retina_red_avg = pd.concat([eid_column, other_columns], axis=1)

    ## Check the number of subjects:
    print("Check the num of subjects for the main ipds")
    main_idps_val = [key for key in dict_final_IDPs.values()]
    for name in main_idps_val:
        count = df_merge_eid_retina_red_avg[name].count()
        print(name, count)

    print(f"Saving file ... {DIR_SAVED_NAME_FILE}")
    df_merge_eid_retina_red_avg.to_csv(DIR_SAVED_NAME_FILE, header=True, index=False)

    # print(df_merge_eid_retina_red_avg['Mean CBF in Putamen (R)'].count())
    # print(df_merge_eid_retina_red_avg['Mean intensity of vessel (L)'].count())
    # print(df_merge_eid_retina_red_avg['Mean intensity of vessel (R)'].count())
    # print(df_merge_eid_retina_red_avg['LA ejection fraction'].count())
