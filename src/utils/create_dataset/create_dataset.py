import pandas as pd
import sys
import logging

from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)
sys.path.append(DIR_UTILS)
from utils.data_information.data_retina_info import list_retina_IDPs, list_retina_IDPs_new, list_retina_homologous, list_retina_homologous_new, list_retina_homologous_red, list_retina_homologous_red_new 


def search_incomplete_column_names(df, list_key):
    """
    Searches for incomplete column names in a DataFrame based on a list of key values.

    Args:
        df (pandas.DataFrame): The DataFrame to search in.
        list_key (list): A list of key values to search for.

    Returns:
        list: A sorted list of column names that match the search criteria.
    """
    filtered_columns = set()

    for datafield in list_key:
        # Filter columns where the number is before the hyphen
        filtered = [col for col in df.columns if col.startswith(f'{datafield}-')]
        filtered_columns.update(filtered)
    return sorted(filtered_columns)

def show_file_and_datafields(DIR_UKBB, list_folders_files, list_keys):
    """
    Displays the file names and the corresponding data fields found in each file.

    Args:
        DIR_UKBB (str): The directory path where the files are located.
        list_folders_files (list): A list of file names to process.
        list_keys (list): A list of data field names to search for.

    Returns:
        tuple: A tuple containing two dictionaries. The first dictionary contains the file names as keys and a list of
        data fields found in each file as values. The second dictionary contains the data field names that were not found
        in any of the files.
    """
    all_columns_found = {}
    columns_not_found = set(list_keys)
    previous_columns = set()

    for file_uk in list_folders_files:
        file_merge_ukbb =  f'{DIR_UKBB}{str(file_uk)}'
        df_file_merge_ukbb = pd.read_csv(file_merge_ukbb, nrows=1)
        list_keys_to_search = set(list_keys) - previous_columns

        list_columns_vascular = search_incomplete_column_names(df_file_merge_ukbb, list_keys_to_search)

        all_columns_found[file_uk] = list_columns_vascular
        columns_found = {col.split('-')[0] for col in list_columns_vascular}
        previous_columns.update(columns_found)
        columns_not_found -= columns_found

    return all_columns_found, columns_not_found

def check_phenotypes_values_count(df):
    """
    Check the values in the specified columns of a DataFrame and compare them to expected counts.
    The values are obtained from the UK Biobank data.

    Args:
        df (pandas.DataFrame): The DataFrame to check.

    Returns:
        None

    Prints a message for each column, indicating whether the actual count of values meets the expected count.

    Example:
        check_phenotypes_values_count(df)
    """
    columns_to_check = {
        '22420-2.0': 39614,
        '22420-3.0': 983,
        '22670-2.0': 50141,
        '22670-3.0': 5186,
        '24103-2.0': 39285,
        '24103-3.0': 915,
        '24372-2.0': 3414,
        '24372-3.0': 1800,
        '24373-2.0': 3414,
        '24373-3.0': 1800
    }

    for col, expected_count in columns_to_check.items():
        actual_count = df[col].count()
        if actual_count >= expected_count:
            logger.info(f"Column {col} has {actual_count} values, which meets the expectation of at least {expected_count}.")
        else:
            logger.warning(f"Column {col} has only {actual_count} values, which is less than the expected count of {expected_count}.")

def merge_and_size(df1, df2, key_used, how_used):
    """
    Merge two DataFrames based on a specified key and merge method, and return the merged DataFrame.

    Args:
        df1 (pandas.DataFrame): The first DataFrame to be merged.
        df2 (pandas.DataFrame): The second DataFrame to be merged.
        key_used (str or list of str): The column(s) used as the key for merging.
        how_used (str): The method used for merging. Can be one of 'inner', 'outer', 'left', or 'right'.

    Returns:
        pandas.DataFrame: The merged DataFrame.

    """
    df_merge = df1.merge(df2, on=key_used, how=how_used)
    logger.info('len(df1), len(df2), len(df_merge)', len(df1), len(df2), len(df_merge))    
    return df_merge

def get_retina_pheno_names(IDPs_retina_used):
    """
    Get the names and identifiers for the retina IDPs based on the specified type.

    Parameters:
    - IDPs_retina_used (str): The type of retina IDPs to use. Possible values are 'multitrait', 'homologous_all', and 'homologous_filtered'.

    Returns:
    - list_retina (list): The list of retina IDPs based on the specified type.
    - list_retina_new (list): The updated list of retina IDPs based on the specified type.
    - retina_eid_identifier (str): The identifier for the retina IDPs.
    - IDPs_retina_name (str): The name of the retina IDPs used.
    """
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

