import itertools
import numpy as np
import pandas as pd
import re
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr, kendalltau, spearmanr
import os, sys
import warnings
import logging

from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)

sys.path.append(DIR_UTILS)
from utils.data_information.data_retina_info import list_retina_IDPs, list_retina_IDPs_new, list_retina_homologous, list_retina_homologous_new, list_retina_homologous_red, list_retina_homologous_red_new
from utils.data_information.data_IDPs_info import list_keys_brain,list_values_brain, list_keys_heart,list_values_heart, list_keys_carotid,list_values_carotid
from utils.data_information.data_homologous_info import dict_homologous, dict_homologous_averag

### Selection datafields:
def custom_sort(column_name):
    """
    Custom sort function for sorting column names.

    Args:
        column_name (str): The column name to be sorted.

    Returns:
        tuple: A tuple containing the sorted column name.

    """
    parts = column_name.split('-')
    if len(parts) == 2 and parts[0].isdigit() and parts[1].replace('.', '', 1).isdigit():
        return (int(parts[0]), float(parts[1]))
    return (column_name,)

def search_incomplete_column_names(df, list_key, list_value):
    """
    Searches for incomplete column names in a DataFrame based on a list of key values.

    Args:
        df (pandas.DataFrame): The DataFrame to search in.
        list_key (list): A list of key values to search for.
        list_value (list): A list of corresponding values for each key.

    Returns:
        list: A sorted list of column names that match the search criteria.
    """
    filtered_columns = set()
    for datafield in list_key:
        # Filter columns where the number is before the hyphen
        filtered = [col for col in df.columns if col.startswith(f'{datafield}-')]
        filtered_columns.update(filtered)
    return sorted(filtered_columns)

def show_file_and_datafields(DIR_UKBB, list_folders_files, list_keys, list_values):
    """
    Displays the file names and corresponding data fields found in the given list of files.

    Parameters:
    - DIR_UKBB (str): The directory path where the files are located.
    - list_folders_files (list): A list of file names to search for.
    - list_keys (list): A list of data field keys to search for.
    - list_values (list): A list of corresponding data field values to search for.

    Returns:
    - all_columns_found (dict): A dictionary containing the file names as keys and the corresponding
                                data fields found as values.
    - columns_not_found (set): A set of data field keys that were not found in any of the files.
    """
    all_columns_found = {}
    columns_not_found = set(list_keys)
    previous_columns = set()

    for file_uk in list_folders_files:
        file_merge_ukbb =  f'{DIR_UKBB}{str(file_uk)}'
        df_file_merge_ukbb = pd.read_csv(file_merge_ukbb, nrows=1)
        list_keys_to_search = set(list_keys) - previous_columns

        list_columns_vascular = search_incomplete_column_names(df_file_merge_ukbb, list_keys_to_search, list_values)

        all_columns_found[file_uk] = list_columns_vascular
        columns_found = {col.split('-')[0] for col in list_columns_vascular}
        previous_columns.update(columns_found)
        columns_not_found -= columns_found

    return all_columns_found, columns_not_found


def print_incomplete_column_names(df, list_key, list_value):
    """
    Prints the names of incomplete columns in a DataFrame.

    Args:
        df (pandas.DataFrame): The DataFrame to check for incomplete columns.
        list_key (list): A list of datafield names.
        list_value (list): A list of corresponding values for the datafields.

    Returns:
        None
    """
    for idx, datafield in enumerate(list_key):
        spike_cols = [col for col in df.columns if datafield in col]
        
        if not spike_cols:
            logger.info(f'Datafield {datafield} is missing, corresponding value: {list_value[idx]}')


def create_counts_df(df):
    """
    Create a DataFrame that stores the counts of non-null values between columns of the input DataFrame.

    Parameters:
        df (pandas.DataFrame): The input DataFrame.

    Returns:
        pandas.DataFrame: A DataFrame with counts of non-null values between columns.
    """
    columns = df.columns
    num_columns = len(columns)

    # Create an empty DataFrame to store the counts
    df_number = pd.DataFrame(index=columns, columns=columns, dtype=float)

    for i, j in itertools.product(range(num_columns), range(num_columns)):
        notna_count = np.sum(df[columns[i]].notna() & df[columns[j]].notna())
        df_number.iloc[i, j] = float(notna_count)

    return df_number #df_number.to_numpy().astype(float)

    
def check_instances(df):
    """
    Check the instances in the given DataFrame.

    Args:
        df (pandas.DataFrame): The DataFrame to check.

    Raises:
        SystemExit: If the instances are different from {'3.0', '2.0'}.

    Returns:
        None
    """
    pattern = re.compile(r'-(\d+\.\d+)')
    versions = set()
    for col in df.columns:
        match = pattern.search(col)
        if match:
            versions.add(match.group(1))
    logger.info('Instances found:', versions)
    if versions != {'3.0', '2.0'}:
        logger.error('The values are different from {3.0, 2.0}. This should not be the case.')
        sys.exit(1)
    else:
        logger.info('Instances are correct. Continuing...')


def combine_instances_23(df):
    """
    Combine instances in the given DataFrame by averaging values for the same parameter and different versions.

    Args:
        df (pandas.DataFrame): The DataFrame containing the data to be processed.

    Returns:
        pandas.DataFrame: The DataFrame with instances combined and averaged.

    """
    versions = set()
    for col in df.columns:
        match = re.search(r'-(\d+\.\d+)', col)
        if match:
            versions.add(match.group(1))

    param_names = {col.split('-')[0] for col in df.columns if '-' in col}
    for param in param_names:
        averaged_rows = 0  # Counter for the number of rows where averaging was performed
        
        for version in versions:
            col_2 = f'{param}-{version}'
            col_3 = f'{param}-{version}'
            
            mask_2 = df[col_2].notna()
            mask_3 = df[col_3].notna()
            
            # Convert columns to numeric if they are not already
            df[col_2] = pd.to_numeric(df[col_2], errors='coerce')
            df[col_3] = pd.to_numeric(df[col_3], errors='coerce')
            
            both_notna = mask_2 & mask_3
            num_both_notna = both_notna.sum()
            
            df[param] = np.where(mask_2, df[col_2], df[col_3])
            df.loc[both_notna, param] = (df[col_2] + df[col_3]) / 2
            
            averaged_rows += num_both_notna
        
        logger.info(f'Number of rows where averaging was performed for {param}: {averaged_rows}')
    
    return df, param_names

def check_combination_higher_indiv(df, var):
    """
    Check if the combined column has at least as many non-null values as the most populated individual column.

    Args:
        df (pandas.DataFrame): The DataFrame containing the data.
        var (str): The name of the variable/column to check.

    Returns:
        None
    """
    notnull_0_0 = df[f'{var}-0.0'].notnull().sum()
    notnull_1_0 = df[f'{var}-1.0'].notnull().sum()
    notnull_2_0 = df[f'{var}-2.0'].notnull().sum()
    notnull_combined = df[var].notnull().sum()
    
    # Make sure that notnull_combined >= max of the individuals
    max_notnull = max(notnull_0_0, notnull_1_0, notnull_2_0)
    if notnull_combined >= max_notnull:
        logger.info(f'Continue: The combined column {var} has at least as many non-null values as the most populated individual column.')
    else:
        logger.error(f'Error: The combined column {var} has fewer non-null values ({notnull_combined}) than the most populated individual column ({max_notnull}).')

def assert_combined_notnull_higher(df, combined_col, individual_cols):
    """
    Asserts that the combined column has a higher number of non-null values than any of the individual columns.

    Args:
        df (pandas.DataFrame): The DataFrame containing the columns.
        combined_col (str): The name of the combined column.
        individual_cols (list): A list of names of individual columns.

    Raises:
        AssertionError: If the combined column has fewer non-null values than the most populated individual column.

    """
    notnull_combined = df[combined_col].notnull().sum()
    notnull_individuals = [df[col].notnull().sum() for col in individual_cols]
    max_notnull_individual = max(notnull_individuals)
    logger.info(f'Checking {combined_col}: {notnull_combined}>={max_notnull_individual}')
    assert notnull_combined >= max_notnull_individual, \
        f'Error: The combined column {combined_col} has fewer non-null values ({notnull_combined}) than the most populated individual column ({max_notnull_individual}).'


def covariants_preprocessing(df):
    # Sex '31-0.0' a '31'
    df.rename(columns={'31-0.0': '31'}, inplace=True)
    # Standing heigh, complete la columna '50' usando '50-0.0', '50-1.0' y '50-2.0' si es necesario
    df['50'] = df['50-2.0'].combine_first(df['50-1.0']).combine_first(df['50-0.0'])
    check_combination_higher_indiv(df, '50')

    # DBP averag '4079-2.0' and '4079-2.1'
    df['4079'] = df[['4079-2.0', '4079-2.1']].mean(axis=1)
    assert_combined_notnull_higher(df, '4079', ['4079-2.0', '4079-2.1'])

    # SBP averag '4080-2.0' and '4080-2.1'
    df['4080'] = df[['4080-2.0', '4080-2.1']].mean(axis=1)
    assert_combined_notnull_higher(df, '4080', ['4080-2.0', '4080-2.1'])

    # Age recruitment 21003
    df['21003'] = df['21003-2.0'].combine_first(df['21003-1.0']).combine_first(df['21003-0.0'])
    check_combination_higher_indiv(df, '21003')

    # Assesment center '54' combining '54-0.0', '54-1.0' y '54-2.0'
    df['54'] = df['54-2.0'].combine_first(df['54-1.0']).combine_first(df['54-0.0'])
    check_combination_higher_indiv(df, '54')

    return df

def adding_extra_cov(df):
    """
    Adds extra covariates to the given DataFrame.

    Args:
        df (pandas.DataFrame): The input DataFrame.

    Returns:
        pandas.DataFrame: The DataFrame with extra covariates added.
    """
    SBP_threshold, DBP_threshold = 140, 90
    
    df['Age_2'] = df['Age']**2
    df['MAP'] = df['DBP']+ 1/3*(df['SBP']-df['DBP'])
    df['Hypertense'] = (df['SBP'] >= SBP_threshold) | (df['DBP'] >= DBP_threshold)
    return df

def check_columns_exist(df, cols):
    """
    Checks if columns exist in the DataFrame.

    Args:
    - df: DataFrame to check columns in.
    - cols: List of column names to check for existence.

    Returns:
    - List of column names that do not exist in the DataFrame.
    """
    return [col for col in cols if col not in df.columns]

def check_and_filter_columns(df, list_retina, list_keys, list_values, eid):
    """
    Checks for missing columns and filters columns in the DataFrame.

    Args:
    - df: DataFrame to filter columns from.
    - list_retina: List of retina column names.
    - list_keys: List of keys column names.
    - list_values: List of values column names.
    - eid: Name of the eid column.

    Returns:
    - Filtered DataFrame with necessary columns.
    """
    try:
        missing_cols = check_columns_exist(df, list_retina + list_keys)
        if missing_cols:
            warnings.warn(f"Warning: Some columns do not exist in the DataFrame: {missing_cols}", UserWarning)

        df_merge_eid_retina_red = df[[eid] + list_retina + list_keys].copy()
    except KeyError as e:
        warnings.warn(f"Warning: Some columns do not exist in the DataFrame: {e}", UserWarning)
        return None

    return df_merge_eid_retina_red

def rename_columns(df, list_keys, list_values, list_retina, list_retina_new):
    """
    Renames columns in the DataFrame.

    Args:
    - df: DataFrame to rename columns in.
    - list_keys: List of keys column names.
    - list_values: List of values column names.
    - list_retina: List of retina column names.
    - list_retina_new: List of new retina column names.

    Returns:
    - DataFrame with renamed columns.
    """
    column_rename_dict = dict(zip(list_keys, list_values))
    column_rename_dict2 = dict(zip(list_retina, list_retina_new))

    return df.rename(columns=column_rename_dict).rename(
        columns=column_rename_dict2
    )

def check_final_columns(df_merge_eid_retina_red, list_retina_new, list_values, eid):
    """
    Checks for missing columns in the final DataFrame.

    Args:
    - df_merge_eid_retina_red: DataFrame to check for missing columns.
    - list_retina_new: List of new retina column names.
    - list_values: List of values column names.
    - eid: Name of the eid column.

    Returns:
    - None
    """
    try:
        missing_cols = check_columns_exist(df_merge_eid_retina_red, list_retina_new + list_values + [eid])
        if missing_cols:
            warnings.warn(f"Warning: Some columns do not exist in the DataFrame: {missing_cols}", UserWarning)
    except KeyError as e:
        warnings.warn(f"Warning: Some columns do not exist in the DataFrame: {e}", UserWarning)



def check_column_counts(df, base_column, comparison_columns, min_count=None):
    base_count = df[base_column].count()
    logger.info(f'{base_column}: {base_count}')

    if min_count is not None and base_count < min_count:
        logger.info(f'Error in {base_column}: count {base_count} is less than {min_count}')

    for comp_col in comparison_columns:
        comp_count = df[comp_col].count()
        logger.info(f'{comp_col}: {comp_count}')
        
        if comp_count < base_count:
            logger.error(f'{comp_col} count {comp_count} is less than {base_column} count {base_count}')

def combine_IMT_cols(df):
    """
    Combines columns containing 'IMT' and 'Min', 'Max', or 'Mean' in their names to calculate the mean IMT values.
    
    Args:
        df (DataFrame): The input DataFrame containing the IMT columns.
        
    Returns:
        DataFrame: The input DataFrame with additional columns for the mean IMT values.
    """
    contains_imt_and_min = [col for col in df.columns if 'IMT' in col and 'Min' in col]
    contains_imt_and_max = [col for col in df.columns if 'IMT' in col and 'Max' in col]
    contains_imt_and_mean = [col for col in df.columns if 'IMT' in col and 'Mean' in col]

    if contains_imt_and_min:
        df['Min carotid IMT'] = df[contains_imt_and_min].mean(axis=1)

    if contains_imt_and_max:
        df['Max carotid IMT'] = df[contains_imt_and_max].mean(axis=1)

    if contains_imt_and_mean:
        df['Mean carotid IMT'] = df[contains_imt_and_mean].mean(axis=1)

    count_contains_max = df[contains_imt_and_max].count().max()
    count_contains_mean = df[contains_imt_and_mean].count().max()
    count_contains_min = df[contains_imt_and_min].count().max()

    count_max = df['Max carotid IMT'].count()
    count_mean = df['Mean carotid IMT'].count()
    count_min = df['Min carotid IMT'].count()

    logger.info('count_max, count_contains_max', count_max, count_contains_max)
    logger.info('count_mean, count_contains_mean', count_mean, count_contains_mean)
    logger.info('count_min, count_contains_min', count_min, count_contains_min)
    return df

def combine_CBF_ATT_cols(df):
    contains_CBF = [col for col in df.columns if 'CBF' in col]
    contains_ATT = [col for col in df.columns if 'ATT' in col]

    if contains_CBF:
        df['Mean CBF'] = df[contains_CBF].mean(axis=1)

    if contains_ATT:
        df['Mean ATT'] = df[contains_ATT].mean(axis=1)

    # Calculate the number of non-NaN values in the previous columns
    count_contains_cbf = df[contains_CBF].count().max()
    count_contains_att = df[contains_ATT].count().max()

    ### print
    # Calculate the number of non-NaN values in the new columns
    count_cbf = df['Mean CBF'].count()
    count_att = df['Mean ATT'].count()

    logger.info('count_cbf, count_contains_cbf', count_cbf, count_contains_cbf)
    logger.info('count_att, count_contains_att', count_att, count_contains_att)

    return df

def analyze_variable_distributions(df):
    '''Analyze variable distributions and differences between time points.

    Args:
        df (pd.DataFrame): The input DataFrame with versioned columns.

    Returns:
        summary_df (pd.DataFrame): Summary DataFrame with non-NaN counts, means, stds, and differences.
    '''
    # Extract unique variables and times
    variables = set()
    times = set()
    for col in df.columns:
        match = re.search(r'(.+)-(\d+\.\d+)', col)
        if match:
            variables.add(match.group(1))
            times.add(float(match.group(2)))

    times = sorted(times)

    summary_data = []

    for var in variables:
        var_data = {}
        var_data['variable'] = var

        # Collect statistics for each time point
        for time in times:
            col_name = f'{var}-{time}'
            if col_name in df.columns:
                non_nan_count = df[col_name].notna().sum()
                mean_val = df[col_name].mean()
                #std_val = df[col_name].std()

                var_data[f'count_time_{time}'] = non_nan_count
                var_data[f'mean_time_{time}'] = mean_val
                #var_data[f'std_time_{time}'] = std_val

        # Check for data at both time points and calculate differences
        col_time_2 = f'{var}-2.0'
        col_time_3 = f'{var}-3.0'

        if col_time_2 in df.columns and col_time_3 in df.columns:
            both_notna_mask = df[col_time_2].notna() & df[col_time_3].notna()
            both_notna_count = both_notna_mask.sum()

            if both_notna_count > 0:
                differences = df.loc[both_notna_mask, col_time_2] - df.loc[both_notna_mask, col_time_3]
                mean_diff = differences.mean()
                #std_diff = differences.std()

                var_data['count_both_times'] = both_notna_count
                var_data['mean_difference'] = mean_diff
                #var_data['std_difference'] = std_diff
                var_data['count_time_3.0_minus_both_times'] = (
                    var_data['count_time_3.0'] - both_notna_count
                )

        # Calculate the total count of non-NaN values for both times
        var_data['total_count'] = var_data.get('count_time_2.0', 0) + var_data.get('count_time_3.0', 0)

        summary_data.append(var_data)

    summary_df = pd.DataFrame(summary_data)
    # Sort the summary DataFrame by count_time_3.0_minus_both_times in descending order
    summary_df = summary_df.sort_values(by='count_time_3.0_minus_both_times', ascending=False)
    return summary_df

def compare_instances_2_3(df_vascular, DIR_SAVED_NAME_STAS):
    threshold_2_3 = 300

    summary_df = analyze_variable_distributions(df_vascular)
    logger.info(summary_df.head(3))
    #save csv
    summary_df.to_csv(DIR_SAVED_NAME_STAS)
    
    max_3_2=int(summary_df['count_time_3.0_minus_both_times'].max())
    logger.info(f'For all the idps, the one with more participants in instance 3.0 has {max_3_2} new ones regarding its instance 2.0')
    if max_3_2<threshold_2_3:
        logger.info(f'Since the number of addings is less than {threshold_2_3}, we just continue with instances 2')

# Filter by specific rows (index) and columns
def filter_dataframe(df, rows, columns):
    '''_summary_

    Args:
        df (_type_): _description_
        rows (_type_): _description_
        columns (_type_): _description_

    Returns:
        _type_: _description_
    '''
    missing_rows = [row for row in rows if row not in df.index]
    missing_columns = [col for col in columns if col not in df.columns]
    
    if missing_rows:
        logger.info(f'Warning: The following rows do not exist in the DataFrame: {missing_rows}')
    if missing_columns:
        logger.info(f'Warning: The following columns do not exist in the DataFrame: {missing_columns}')
    
    # Filtrar solo las etiquetas que existen
    rows = [row for row in rows if row in df.index]
    columns = [col for col in columns if col in df.columns]
    
    return df.loc[rows, columns]


def get_retina_names(IDPs_retina_used):
    '''_summary_

    Args:
        IDPs_retina_used (_type_): _description_

    Returns:
        _type_: _description_
    '''
    ### RETINA IDP
    if IDPs_retina_used == 'multitrait':
        list_retina = list_retina_IDPs
        list_retina_new = list_retina_IDPs_new
        retina_eid_identifier = '0'
        IDPs_retina_name = IDPs_retina_used
        
    elif IDPs_retina_used == 'homologous_all':
        list_retina = list_retina_homologous
        list_retina_new = list_retina_homologous_new
        retina_eid_identifier = 'eid'
        IDPs_retina_name = IDPs_retina_used

    elif IDPs_retina_used=='homologous_filtered':
        logger.info('Filter made, only homologous filtered') #based on the csv
        list_retina = list_retina_homologous_red
        list_retina_new = list_retina_homologous_red_new
        retina_eid_identifier = 'eid'
        IDPs_retina_name = 'homologous_all'
    
    return list_retina, list_retina_new, retina_eid_identifier, IDPs_retina_name


def get_IDPs_names(IDPs_vascular_used, combine=False):
    '''_summary_

    Args:
        IDPs_vascular_used (_type_): _description_
        combine_IMT (_type_): _description_

    Returns:
        _type_: _description_
    '''
    aveg_keys = list(dict_homologous_averag.keys())
    aveg_values = list(dict_homologous_averag.values())
    assert len(aveg_keys) == len(aveg_values), "The lengths of aveg_keys and aveg_values are not equal."

    ### VASCULAR IDP
    if (IDPs_vascular_used=='all') and (combine!=False):
        list_keys = list_keys_brain + list_keys_heart + list_keys_carotid + aveg_keys
        list_values = list_values_brain + list_values_heart + list_values_carotid + aveg_values
        ## TODO, they are repeted righ not!
        list_keys=list(set(list_keys))
        list_values=list(set(list_values))
        logger.info('All UKBB IDPs, also homologous combined')

    elif IDPs_vascular_used == 'all':
        list_keys = list_keys_brain + list_keys_heart + list_keys_carotid 
        list_values = list_values_brain + list_values_heart + list_values_carotid
        logger.info('All UKBB IDPs, NO homologous combined')

    elif IDPs_vascular_used=='homologous_all':
        logger.info('Filter made, only homologous') #based on the csv
        if combine:
            list_keys = aveg_keys
            list_values = aveg_values
        else:
            list_keys = list(dict_homologous.keys())
            list_values = list(dict_homologous.values())


    assert len(list_keys) == len(list_values), "The lengths of list_keys and list_values are not equal."

    return list_keys, list_values


def get_phenotypes_names(IDPs_retina_used, IDPs_vascular_used, combine):
    '''_summary_

    Args:
        IDPs_retina_used (_type_): _description_
        IDPs_vascular_used (_type_): _description_
        combine (_type_): _description_

    Returns:
        _type_: _description_
    '''
    list_retina, list_retina_new, retina_eid_identifier, IDPs_retina_name = get_retina_names(IDPs_retina_used)
    list_keys, list_values = get_IDPs_names(IDPs_vascular_used, combine)

    return list_retina, list_retina_new, retina_eid_identifier, IDPs_retina_name, list_keys, list_values

def complete_column_instaces(df, list_new_columns):
    for new_column in list_new_columns:
        start_new_column = f'{new_column}-'
        df[new_column] = None

        # Llenar la nueva columna 'X' con los valores de las columnas que siguen el patrón '54-*'
        for col in df.columns:
            if col.startswith(start_new_column):
                df[new_column] = df[new_column].combine_first(df[col])
        
        # Verificar que df['X'] tenga al menos tantos valores no nulos como la columna 'X-' que más tenga
        assert df[new_column].count() >= max(df[col].count() for col in df.columns if col.startswith(start_new_column))

        
        # Eliminar las columnas originales que siguen ese patrón si es necesario
        df = df.drop(columns=[col for col in df.columns if col.startswith(start_new_column)])
    return df

#def kendall_pval(x,y):
#    return kendalltau(x,y)[1]
#def spearmanr_pval(x,y):
#    return spearmanr(x,y)[1]

def pearsonr_pval(x,y):
    return pearsonr(x,y)[1]

def check_vals(df_used, dict_final_IDPs):
    logger.info('Check the num of subjects for the main ipds')
    main_idps = [val for val in dict_final_IDPs.values()if val != 'Min carotid IMT']
    logger.info(main_idps)
    for name in main_idps:
        count = df_used[name].count()
        logger.info(name, count)

def check_keys(df_used, dict_final_IDPs):
    logger.info('Check the num of subjects for the main ipds')
    main_idps = [key for key in dict_final_IDPs.keys() if key != 'Min carotid IMT']
    for name in main_idps:
        count = df_used[name].count()
        logger.info(name, count)

def replace_outliers_with_nan(df_, columns, NUM_STD, dir_to_save):
    stats = []

    for col in columns:
        mean = df_[col].mean()
        std = df_[col].std()
        upper_limit = mean + NUM_STD * std
        lower_limit = mean - NUM_STD * std
        
        non_nan_before = df_[col].notna().sum()
        df_.loc[(df_[col] > upper_limit) | (df_[col] < lower_limit), col] = np.nan
        
        non_nan_after = df_[col].notna().sum()
        
        stats.append([col, mean, std, non_nan_before, non_nan_after])
    
    stats_df = pd.DataFrame(stats, columns=['Column', 'Mean', 'STD', 'Non_NaN_Before', 'Non_NaN_After'])
    stats_df.to_csv(f"{dir_to_save}remove_outliers_column_stats_std_{NUM_STD}.csv", index=False)
    
    return df_

def replace_outliers_IQR_with_nan(df, columns, save_dir):
    """
    Replaces outliers in specified columns of the DataFrame with NaN using the IQR method.
    
    Parameters:
    df (pd.DataFrame): The input DataFrame.
    columns (list): List of column names to check for outliers.
    save_dir (str): Directory to save the DataFrame with NaNs.
    
    Returns:
    pd.DataFrame: DataFrame with outliers replaced by NaN.
    """
    
    df_out = df.copy()
    
    for col in columns:
        if col in df_out.columns:
            Q1 = df_out[col].quantile(0.25)
            Q3 = df_out[col].quantile(0.75)
            IQR = Q3 - Q1
            lower_bound = Q1 - 1.5 * IQR
            upper_bound = Q3 + 1.5 * IQR
            df_out[col] = df_out[col].apply(lambda x: x if lower_bound <= x <= upper_bound else np.nan)
        else:
            print(f"Column {col} not found in DataFrame")

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    file_path = os.path.join(save_dir, "df_with_nans.csv")
    df_out.to_csv(file_path, index=False)
    print(f"DataFrame saved to {file_path}")
    
    return df_out

def check_column_counts(df, min_count=3000):
    """
    Check non-null counts of DataFrame columns and raise an error if any column has fewer than min_count non-null values.

    Parameters:
    df (pd.DataFrame): DataFrame to check.
    min_count (int): Minimum number of non-null values required for each column.

    Returns:
    list: List of tuples containing column names and their non-null counts.
    """
    column_counts = []
    
    for col in df.columns:
        count = df[col].count()
        column_counts.append((col, count))
        
        if count < min_count:
            raise ValueError(
                f"Error: The column {col} has fewer than {min_count} non-null values ({count})"
            )
    
    return column_counts


def check_column_thresholds(df, column_thresholds):
    """
    Check if columns in a DataFrame meet specified non-null value thresholds.

    Parameters:
    df (pd.DataFrame): DataFrame to check.
    column_thresholds (dict): Dictionary with column names as keys and minimum count thresholds as values.

    Raises:
    ValueError: If any column has fewer non-null values than the specified threshold.
    """
    for column, threshold in column_thresholds.items():
        count = df[column].count()
        
        if count < threshold:
            raise ValueError(f"{column} error! {count} is < {threshold}")
        else:
            print(f"{column} {count} is >= {threshold}")


def get_baselines(df, list_main_IDPs):
    """
    Calculate baselines for a given DataFrame and a list of main IDPs.

    Parameters:
    - df (pandas.DataFrame): The input DataFrame containing the data.
    - list_main_IDPs (list): The list of main IDPs (Independent Data Points).

    Returns:
    - df_baselines (pandas.DataFrame): The DataFrame containing the calculated baselines.

    """
    baselines = []
    cols_used = list_main_IDPs + ['Age', 'Sex']
    df_subset = df[cols_used]  
    for idp in list_main_IDPs:
        df_temp = df_subset[[idp, 'Age', 'Sex']].dropna()
        non_null = df_temp.shape[0]
        age_mean = df_temp['Age'].mean()
        age_std = df_temp['Age'].std()
        sex_counts = df_temp['Sex'].value_counts().to_dict()  # Convert to dictionary
        
        # Append a dictionary with appropriate keys
        baselines.append({
            'IDP': idp,
            'non_null': non_null,
            'Age_Mean': age_mean,
            'Age_Std': age_std,
            'Sex_Counts': sex_counts
        })
    
    df_baselines = pd.DataFrame(baselines)
    return df_baselines