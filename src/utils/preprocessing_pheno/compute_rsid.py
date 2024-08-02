import pandas as pd
from sklearn.linear_model import LinearRegression
import sys
import logging

from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)

sys.path.append(DIR_UTILS)

from utils.data_information.data_covariants_info import covariants_GPC, covariants_common
from utils.create_dataset.create_dataset import merge_and_size

def get_cov_names(covariants):
    """
    Returns a list of column names based on the given covariants.

    Args:
        covariants (str): The type of covariants to consider.

    Returns:
        list: A list of column names.

    Raises:
        None

    Examples:
        >>> get_cov_names('age_age2_sex')
        ['Age', 'Sex', 'Age_2']

        >>> get_cov_names('age_age2_sex_PCs')
        ['Age', 'Sex', 'Age_2', 'PC1', 'PC2', 'PC3']

        >>> get_cov_names('all_BP')
        ['Common1', 'Common2', 'PC1', 'PC2', 'Age_2', 'SBP', 'DBP']

    """
    list_PCs = list(covariants_GPC.values())
    if covariants == 'age_age2_sex':
        list_cov_columns = ['Age', 'Sex', 'Age_2']

    elif covariants == 'age_age2_sex_PCs':
        list_cov_columns = ['Age', 'Sex', 'Age_2'] + list_PCs

    elif covariants in ['BP', 'Hypertense', 'MAP']:
        list_new_columns_first = list(covariants_common.values())
        if covariants == 'BP':
            list_BP = ['SBP','DBP']
        elif covariants == 'all_MAP':
            list_BP = ['MAP']
        elif covariants == 'Hypertense':
            list_BP = ['Hypertense']
        list_age_sex_comp = ['Age_2'] #, '21003_31']
        list_cov_columns = list_new_columns_first + list_PCs + list_age_sex_comp + list_BP
    else:
        list_new_columns_first = list(covariants_common.values()) #['21003','31','54','50'] #,'21002','21001']
        #logger.info('sex*age')
        #logger.info(df['31'].value_counts)
        #df['21003_31'] = df['21003']*df['31']
        #logger.info(df['21003_31'].value_counts)
        list_age_sex_comp = ['Age_2'] #, '21003_31']
        list_cov_columns = list_new_columns_first + list_PCs + list_age_sex_comp

    logger.info(f'Selected cov columns {list_cov_columns}')
    return list_cov_columns


def df_IDPs_cov_substracted(df_merge, common_covariates, common_phenotypes):
    """
    Compute the corrected phenotype values by subtracting the effect of common covariates.

    Args:
        df_merge (pandas.DataFrame): The input DataFrame containing the merged data.
        common_covariates (str or list): The common covariates to be used in the regression model.
            If a string is provided, it will be converted to a list.
        common_phenotypes (list): The list of common phenotypes to be corrected.

    Returns:
        pandas.DataFrame: The DataFrame with the corrected phenotype values.

    Raises:
        None

    """

    # Create an empty DataFrame to store the corrected phenotype values
    df_original = df_merge.copy()
    # Remove the columns in common_phenotypes from df_original
    df_original.drop(columns=common_phenotypes, inplace=True)

    df = pd.DataFrame(index=df_merge.index)
    df['eid'] = df_merge['eid']

    # Iterate through each common phenotype
    for phenotype in common_phenotypes:
        logger.info(phenotype)
        if isinstance(common_covariates, str): # for BP
            common_covariates = [common_covariates]
        # Create a temporary DataFrame with the columns of interest for the current pair
        temp_df_aux = df_merge[common_covariates + [phenotype]]

        # Drop rows with NaN values in any of the columns
        temp_df = temp_df_aux.dropna()

        # Convert all columns to numeric type with errors='coerce'
        temp_df = temp_df.apply(pd.to_numeric, errors='coerce')

        # Split the data into predictors (X) and target (y) for the current pair
        X_pair = temp_df[common_covariates]
        y_pair = temp_df[phenotype]

        # Check if the '54' variable (if it exist) is categorical
        if 'UK Biobank assessment centre' in common_covariates and temp_df['UK Biobank assessment centre'].dtype == 'O':
            # If it's categorical, encode it as dummy variables
            X_pair = pd.get_dummies(X_pair, columns=['54'], drop_first=True)

        # Create a regression model
        regressor = LinearRegression()

        # Fit the model to the current pair
        regressor.fit(X_pair, y_pair)

        # Calculate the corrected phenotype values
        corrected_phenotype = y_pair - regressor.predict(X_pair)

        # Store the corrected phenotype values in the df DataFrame
        df[phenotype] = corrected_phenotype

        # Check if the length of the corrected phenotype column matches the original length
        if (temp_df[phenotype].notnull().sum()) != (df[phenotype].notnull().sum()):
            logger.info(f'The number of not NaN values in the corrected phenotype column {phenotype} does not match the original number.', (temp_df[phenotype].notnull().sum()), (df[phenotype].notnull().sum()))

        # Check if all values in the current phenotype column are NaN
        if df[phenotype].isnull().all():
            logger.info(f'The phenotype column {phenotype} only contains NaN values.')

    return merge_and_size(df_original, df, 'eid', how_used='left')


def normalize_df(df_merge):
    """
    Normalize the given DataFrame by subtracting the mean and dividing by the standard deviation.

    Parameters:
    - df_merge (pandas.DataFrame): The DataFrame to be normalized.

    Returns:
    - pandas.DataFrame: The normalized DataFrame.
    """
    return df_merge.apply(lambda x: (x - x.mean()) / x.std() if x.notna().any() else x, axis=0)


def save_cov_files(DIR_NEW_NAME_FILE, df, df_merge_eid_not_zcored, list_IDPs):
    """
    Save the covariance files.

    Args:
        DIR_NEW_NAME_FILE (str): The directory and filename to save the file.
        df (pandas.DataFrame): The DataFrame containing the data.
        df_merge_eid_not_zcored (pandas.DataFrame): The DataFrame containing the merged data.
        list_IDPs (list): The list of IDPs.

    Returns:
        None
    """
    ### df has no index, so you have to bring it back from df_merge_eid_not_zcored
    df_save = df.merge(df_merge_eid_not_zcored[['eid']], left_index=True, right_index=True, how='right')
    df_save[['eid'] + list_IDPs]
    logger.info('the two first should be the same: ', len(df_save), len(df_merge_eid_not_zcored), len(df))
    if len(df_save)!=len(df_merge_eid_not_zcored):
        logger.info('ERROR!')
    else:
        logger.info('Ok, saving the file', DIR_NEW_NAME_FILE)
        df_save.to_csv(DIR_NEW_NAME_FILE)


def get_files_names(IDPs_retina_used, IDPs_vascular_used, z_score_applied, combine_IMT_heart_brain, covariants, TYPE_BP, DATE_USED, FILTER_OUTLIERS):
    """
    Get the names of the files based on the input parameters.

    Parameters:
    - IDPs_retina_used (str): The IDPs used for the retina.
    - IDPs_vascular_used (str): The IDPs used for the vascular system.
    - z_score_applied (str): The type of z-score applied.
    - combine_IMT_heart_brain (str): The method used to combine IMT, heart, and brain.
    - covariants (str): The covariants used.
    - TYPE_BP (bool or str): The type of blood pressure. If False, all blood pressure types are used. If 'all_Hyper', only hypertense blood pressure is used.
    - DATE_USED (str): The date used.

    Returns:
    - name_file_values_count (str): The name of the file containing the values count.
    - name_file_zscored (str): The name of the file containing the z-scored values.
    """
    if (TYPE_BP == False) or (TYPE_BP == 'False'):
        common_part = f'{IDPs_retina_used}_{IDPs_vascular_used}_zscored_{z_score_applied}_combine_{combine_IMT_heart_brain}_cov_used_{covariants}_outliers_{FILTER_OUTLIERS}.csv'
    elif TYPE_BP == 'Hypertense': #['BP', 'Hypertense', 'MAP']
        common_part = f'{IDPs_retina_used}_{IDPs_vascular_used}_zscored_{z_score_applied}_combine_{combine_IMT_heart_brain}_cov_used_Hypertense_outliers_{FILTER_OUTLIERS}.csv'
    else:
        raise ValueError(f'Error, TYPE_BP={TYPE_BP} is not defined!')

    name_file_values_count = f'{DATE_USED}_N_subjects_{common_part}'
    name_file_zscored = f'{DATE_USED}_{common_part}'
    return name_file_values_count, name_file_zscored
