import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import numpy as np
import os, glob, sys
import re
import logging

from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS, DIR_LDSR, DIR_OUTPUT, DIR_OTHER_GENES_PATH, DIR_GENES_RET, P_VAL_GENES, P_VAL_PATHS 

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)
sys.path.append(DIR_UTILS)

from utils.plotting.generate_plots import get_size_val
from utils.plotting.plotting_settings import configure_global_params

configure_global_params()

def plot_top_elements_frequency(series, top_n=10):
    """
    Plot the frequency of the top N elements in a series.

    Parameters:
    - series: pandas Series or DataFrame
        The series or dataframe containing the elements to be counted.
    - top_n: int, optional
        The number of top elements to be plotted. Default is 10.

    Returns:
    None
    """
    # Convert the DataFrame to a series for easy value counting
    flattened_series = series.values.flatten()

    # Initialize an empty dictionary for element counting
    element_count = {}

    # Iterate over each set in the series to count unique elements
    for item_set in flattened_series:
        for element in item_set:
            element_count[element] = element_count.get(element, 0) + 1

    # Convert the dictionary into a series for counting
    flattened_series = pd.Series(element_count)

    # Take the top N elements (you can change this number as needed)
    top_values = flattened_series.sort_values(ascending=False).head(top_n)

    # Create the bar chart
    plt.figure(figsize=(10, 6))
    top_values.plot(kind='bar', color='skyblue')
    plt.title('Frequency of Top Elements')
    plt.xlabel('Elements')
    plt.ylabel('Frequency')
    plt.xticks(rotation=45, ha='right')  # Rotate x-axis names for better readability
    plt.tight_layout()
    plt.show()

def max_min_num_genes(df):
    """
    Calculate the minimum and maximum values of a DataFrame and return them along with the normalization method used.

    Args:
        df (pandas.DataFrame): The DataFrame containing the data.

    Returns:
        tuple: A tuple containing the minimum value, maximum value, and normalization method used.

    """
    vmin_val = df.min().min()
    vmax_val = df.max().max()
    norm_used = LogNorm()
    #norm_used = LogNorm() if vmax_val > 160 else LogNorm()
    return vmin_val, vmax_val, norm_used

def mask_plot_num_shared_genes(df, cmap_used, ax, label_cbar=None):
    """
    Plot a heatmap of the correlation matrix, showing only the lower triangular part.

    Args:
        df (pandas.DataFrame): The correlation matrix.
        cmap_used (str or colormap): The colormap to use for the heatmap.
        ax (matplotlib.axes.Axes): The axes on which to draw the heatmap.
        label_cbar (str, optional): The label for the colorbar. Defaults to None.

    Returns:
        matplotlib.axes.Axes: The axes on which the heatmap is drawn.
    """
    mask = np.tril(np.ones_like(df.corr(), dtype=bool))
    df_T = df.where(mask)

    vmin_val, vmax_val, norm_used = max_min_num_genes(df_T)

    return sns.heatmap(df_T, 
            annot=True, 
            fmt=".0f", 
            cbar=True, 
            annot_kws={'weight': 'normal'}, 
            vmin=int(vmin_val),
            vmax=int(vmax_val),
            cmap=cmap_used,
            alpha=1.0, 
            cbar_kws={'label': label_cbar}, 
            norm=norm_used,
            mask=~mask,  # Apply the mask to show only the lower triangular part
            ax=ax)

def no_mask_plot_num_shared_genes(df_T, cmap_used, ax, label_cbar=None):
    """
    Plot the number of shared genes using a heatmap without a mask.

    Parameters:
    - df_T (DataFrame): The input DataFrame containing the data to be plotted.
    - cmap_used (str or Colormap): The colormap to be used for the heatmap.
    - ax (Axes): The matplotlib Axes object on which the heatmap will be plotted.
    - label_cbar (str, optional): The label for the colorbar. Defaults to None.

    Returns:
    - ax (Axes): The matplotlib Axes object with the heatmap plotted.

    """
    vmin_val, vmax_val, norm_used = max_min_num_genes(df_T)
    return sns.heatmap(df_T, 
            annot=True, 
            fmt=".0f", 
            cbar=True, 
            annot_kws={'weight': 'normal'}, 
            vmin=int(vmin_val), 
            vmax=int(vmax_val),
            cmap=cmap_used,
            alpha=1.0, 
            cbar_kws={'label': label_cbar}, 
            norm=norm_used,
            ax=ax)

def doble_num_shared_genes(df_corr1, df_corr2, figsize_val, width_ratios_val, title_fig=False, only_half1=False, only_half2=False, cmap_used1='seismic', cmap_used2='RdBu_r'):
    """
    Plot the number of shared genes between two correlation matrices.

    Args:
        df_corr1 (DataFrame): The first correlation matrix.
        df_corr2 (DataFrame): The second correlation matrix.
        figsize_val (tuple): The size of the figure (width, height).
        width_ratios_val (list): The width ratios of the subplots.
        title_fig (str, optional): The title of the figure. Defaults to False.
        only_half1 (bool, optional): Whether to plot only the upper half of the first correlation matrix. Defaults to False.
        only_half2 (bool, optional): Whether to plot only the upper half of the second correlation matrix. Defaults to False.
        cmap_used1 (str, optional): The colormap used for the first correlation matrix. Defaults to 'seismic'.
        cmap_used2 (str, optional): The colormap used for the second correlation matrix. Defaults to 'RdBu_r'.

    Returns:
        None
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=figsize_val, gridspec_kw={'width_ratios': width_ratios_val})

    if only_half1:
        fig1 = mask_plot_num_shared_genes(df_corr1, cmap_used1, ax1)
    
    elif only_half1==False:
        fig1 = no_mask_plot_num_shared_genes(df_corr1, cmap_used1, ax1)

    fig1.set_xticklabels(ax1.get_xticklabels(), rotation = 45, ha='right')
    
    if only_half2:
        fig2 = mask_plot_num_shared_genes(df_corr2, cmap_used2, ax2)

    elif only_half2==False:
        fig2 = no_mask_plot_num_shared_genes(df_corr2, cmap_used2, ax2)

    fig2.set_xticklabels(ax2.get_xticklabels(), rotation = 45, ha='right')

    # Adjust layout for better spacing
    plt.show()

    if title_fig!=False:
        fig.savefig(f'{DIR_OUTPUT}{title_fig}', dpi=300, format='jpg', bbox_inches='tight', pad_inches=0.1)

def plot_num_shared_genes(df_T, val_conversion_factor=2, val_conversion_factor_extra=1, name_fig=False, only_half=False):
    """
    Plots the number of shared genes.

    Args:
        df_T (DataFrame): The input DataFrame containing the data.
        val_conversion_factor (int, optional): The conversion factor for adjusting the figure size based on the length of x-axis labels. Defaults to 2.
        val_conversion_factor_extra (int, optional): The extra conversion factor for adjusting the figure size. Defaults to 1.
        name_fig (str, optional): The name of the figure file to be saved. Defaults to False.
        only_half (bool, optional): Flag indicating whether to plot only half of the data. Defaults to False.

    Returns:
        Figure: The generated figure object.
    """
    # Adjust the figure size based on the length of x-axis labels
    num_columns = get_size_val(df_T, val_conversion_factor)
    num_rows = get_size_val(df_T, val_conversion_factor, False)
    fig, ax = plt.subplots(figsize=(num_columns/val_conversion_factor_extra, num_rows))
    cmap_used= 'Purples'
    
    if only_half:
        fig1 = mask_plot_num_shared_genes(df_T, cmap_used, ax, label_cbar=None)

    else:
        fig1 = no_mask_plot_num_shared_genes(df_T, cmap_used, ax, label_cbar=None)

    plt.xticks(rotation=45, ha='right')
        
    if name_fig!=False:
        fig.savefig(f'{DIR_OUTPUT}{name_fig}', dpi=300, format='jpg', bbox_inches='tight', pad_inches=0.1)
    
    return fig


def calculate_shared_elements_counts(df):
    """
    Calculate the counts of shared elements between pairs of 'name_pheno' entries in a DataFrame.

    Parameters:
    - df (pandas.DataFrame): The input DataFrame containing 'name_pheno' and 'names' columns.

    Returns:
    - df_shared_elements (pandas.DataFrame): A square DataFrame containing the counts of shared elements between pairs of 'name_pheno' entries.
    - df_shared_elements_names (pandas.DataFrame): A square DataFrame containing the shared elements between pairs of 'name_pheno' entries.

    Example usage:
    df = pd.DataFrame({'name_pheno': ['A', 'A', 'B', 'B'], 'names': [['John', 'Alice'], ['Alice', 'Bob'], ['Bob', 'Charlie'], ['Charlie', 'David']]})
    shared_counts, shared_counts_names = calculate_shared_elements_counts(df)
    """

    # Initialize a dictionary to store the counts of shared elements between pairs
    shared_counts = {}
    shared_counts_names = {}

    # Create a list of unique 'name_pheno' entries
    name_phenos = df['name_pheno'].unique()

    # Iterate through each 'name_pheno' pair for comparison
    for name_pheno1 in name_phenos:
        names1 = set(df[df['name_pheno'] == name_pheno1]['names'].iloc[0])

        for name_pheno2 in name_phenos:
            names2 = set(df[df['name_pheno'] == name_pheno2]['names'].iloc[0])

            # Calculate shared elements and store counts in the dictionary
            pair_key = f"{name_pheno1} with {name_pheno2}"
            shared_counts[pair_key] = len(names1.intersection(names2))
            shared_counts_names[pair_key] = names1.intersection(names2)

    # Create a square DataFrame from the shared_counts dictionary
    df_shared_elements = pd.DataFrame(index=name_phenos, columns=name_phenos)
    df_shared_elements_names = pd.DataFrame(index=name_phenos, columns=name_phenos)
    for index in name_phenos:
        for column in name_phenos:
            df_shared_elements.at[index, column] = shared_counts[f"{index} with {column}"]
            df_shared_elements_names.at[index, column] = shared_counts_names[f"{index} with {column}"]

    return df_shared_elements, df_shared_elements_names


def rename_col_index(df, l_diseases_old, l_diseases_new):
    """
    Renames the columns and index of a DataFrame based on the provided lists of old and new names.

    Args:
        df (pandas.DataFrame): The DataFrame to be modified.
        l_diseases_old (list): A list of old column and index names.
        l_diseases_new (list): A list of new column and index names.

    Returns:
        pandas.DataFrame: The modified DataFrame with renamed columns and index.
    """
    # Selecting only the columns and index that match the 'l_diseases_new' list
    selected_columns = [col for col in df.columns if col in l_diseases_old]
    selected_index = [idx for idx in df.index if idx in l_diseases_old]

    # Filtering the DataFrame to retain only the selected columns and index
    df_filtered = df.loc[selected_index, selected_columns]

    # Renaming index and columns
    df_filtered.rename(index=dict(zip(l_diseases_old, l_diseases_new)), inplace=True)
    df_filtered.rename(columns=dict(zip(l_diseases_old, l_diseases_new)), inplace=True)

    return df_filtered

def reorder_index_columns(df, list_n):
    """
    Reorders the index and columns of a DataFrame according to the provided list.

    Parameters:
    df (pandas.DataFrame): The DataFrame to be reordered.
    list_n (list): The list specifying the desired order of index and columns.

    Returns:
    pandas.DataFrame: The reordered DataFrame.
    """
    df = df.reindex(index=list_n, columns=list_n)
    return df

def compute_intersections_csv(p_value_min, filenames, input_dir, save_results, csv_name_all, csv_name_count, csv_name_diagonal, csv_name, csv_genes_name):
    """
    Compute intersections and save results in CSV files.

    Args:
        p_value_min (float): The minimum p-value threshold.
        filenames (list): List of filenames.
        input_dir (str): The directory where the input files are located.
        save_results (str): The directory where the results will be saved.
        csv_name_all (str): The name of the CSV file to save all significant genes.
        csv_name_count (str): The name of the CSV file to save the count of significant genes per phenotype.
        csv_name_diagonal (str): The name of the CSV file to save the number of significant genes per phenotype.
        csv_name (str): The name of the CSV file to save the intersection lengths.
        csv_genes_name (str): The name of the CSV file to save the names of genes in the intersection.

    Returns:
        tuple: A tuple containing the following DataFrames:
            - df_count: The count of significant genes per phenotype.
            - df_guardar_final: The number of significant genes per phenotype.
            - df_save_shapes: The intersection lengths.
            - df_save_intersections: The names of genes in the intersection.
    """
    l_aux = []

    for file in filenames:
        df = pd.read_csv(input_dir+file+'__gene_scores', delimiter='\t', names=['gen', 'p'])
        df['file_col'] = file
        l_aux.append(df)

    df_concat = pd.concat(l_aux)

    df_concat['-log10(p)'] = -np.log10(df_concat['p'])
    y = df_concat[df_concat['-log10(p)'] >= p_value_min]
    df_significant = y.sort_values('-log10(p)', ascending=False)

    df_significant.to_csv(save_results + csv_name_all + '.csv')
    df_count = df_significant['gen'].value_counts().to_frame()
    df_count['ratio_N_pheno'] = df_significant['gen'].value_counts().to_frame() / len(filenames)
    df_count.to_csv(save_results + csv_name_count + '.csv')

    df_guardar = pd.DataFrame(df_significant.groupby(by=['file_col'])['gen'].apply(list))
    df_guardar2 = pd.DataFrame(df_significant.groupby(by=['file_col'])['gen'].count())
    df_guardar_final = df_guardar.merge(df_guardar2, how='inner', on='file_col')
    df_guardar_final.to_csv(save_results + csv_name_diagonal + '.csv')
    df_save_shapes = pd.DataFrame([])
    df_save_intersections = pd.DataFrame([])

    i = 0
    for file in filenames:
        i = i + 1
        genes = df_significant[df_significant['file_col'] == file]['gen']
        genes = genes.to_list()
        l_aux2 = []
        l_aux3 = []

        for j in range(len(filenames)):
            other_file = filenames[j]
            df_intersection = df_significant[(df_significant['file_col'] == other_file) & (df_significant['gen'].isin(genes))]
            save_shapes = df_intersection.shape[0]
            l_aux3.append(save_shapes)

            l_aux2.append(df_intersection['gen'].to_list())

        df = pd.DataFrame({file: l_aux3})
        df_save_shapes = pd.concat([df_save_shapes, df], axis=1)

        df2 = pd.DataFrame({file: l_aux2})
        df_save_intersections = pd.concat([df_save_intersections, df2], axis=1)

    df_save_shapes = df_save_shapes.set_axis(df_save_shapes.columns, axis='index')
    df_save_shapes.to_csv(save_results + csv_name + '.csv')

    df_save_intersections = df_save_intersections.set_axis(df_save_intersections.columns, axis='index')
    df_save_intersections.to_csv(save_results + csv_genes_name + '.csv')

    return df_count, df_guardar_final, df_save_shapes, df_save_intersections

def values_above_threshold(df, threshold=1):
    """
    Find and display values above a given threshold in a DataFrame.

    Parameters:
    - df (pandas.DataFrame): The DataFrame to search for values above the threshold.
    - threshold (float or int, optional): The threshold value. Defaults to 1.

    Returns:
    None
    """
    # Find and display values above the threshold
    above_threshold = df[df > threshold].stack()
    
    #logger.info(f"Values above {threshold} :\n {above_threshold}")
    
    for (row, col), value in above_threshold.iteritems():
        logger.info(f"Index: {row}, Column: {col}, Value: {value}")

def plot_limits(df):
    """
    Calculate the maximum value in the DataFrame and its corresponding index and column.
    Print information about vmin, vmax, and the maximum value's index and column.
    If the maximum value is greater than 1, set vmax to 1 and issue a warning.

    Args:
        df (pandas.DataFrame): The input DataFrame.

    Returns:
        float: The maximum value (vmax) to be used for plotting.
    """
    vmax = abs(df).max().max()
    #max_index = abs(df).stack().idxmax()
    values_above_threshold(df, threshold=1)
    if vmax > 1:
        logger.warning('For this image the cbar has been force to be into -1 and 1')
        vmax = 1
    return vmax

def mask_plot_corr_pval(df_corr, df_pval, cmap_used, ax, label_cbar=None):
    """
    Plot a masked correlation heatmap with p-values.

    Parameters:
    - df_corr (pandas.DataFrame): The correlation matrix.
    - df_pval (pandas.DataFrame): The p-value matrix.
    - cmap_used (str or colormap): The colormap to use for the heatmap.
    - ax (matplotlib.axes.Axes): The axes on which to draw the heatmap.
    - label_cbar (str, optional): The label for the colorbar. Defaults to None.

    Returns:
    - seaborn.heatmap: The heatmap plot.

    """
    mask = np.tril(np.ones_like(df_corr.corr(), dtype=bool))
    df_corr = df_corr.where(mask)
    vmax_val = plot_limits(df_corr)
    return sns.heatmap(df_corr,  
                annot=df_pval, 
                cbar=True,
                fmt="", 
                annot_kws={'weight': 'normal'}, 
                vmin=-vmax_val, 
                vmax=vmax_val, 
                cmap=cmap_used, 
                alpha=1.0, 
                cbar_kws={'label': label_cbar},
                mask=~mask,
                ax=ax)

def no_mask_plot_corr_pval(df_corr, df_pval, cmap_used, ax, label_cbar=None):
    """
    Plot correlation matrix with p-values.

    Parameters:
    - df_corr (pandas.DataFrame): The correlation matrix.
    - df_pval (pandas.DataFrame): The p-value matrix.
    - cmap_used (str or colormap): The colormap to use for the heatmap.
    - ax (matplotlib.axes.Axes): The axes to plot the heatmap on.
    - label_cbar (str, optional): The label for the colorbar. Defaults to None.

    Returns:
    - seaborn.heatmap: The heatmap plot.

    """
    vmax_val = plot_limits(df_corr)
    return sns.heatmap(df_corr,  
                    annot=df_pval, #.T, 
                    cbar=True, #If not False
                    fmt="", 
                    annot_kws={'weight': 'normal'}, 
                    vmin=-round(vmax_val,2), 
                    vmax=round(vmax_val,2), 
                    cmap=cmap_used,
                    alpha=1.0, 
                    cbar_kws={'label': label_cbar},
                    ax=ax)

def doble_plot_corr_pval(df_corr1, df_pval1, df_corr2, df_pval2, figsize_val, width_ratios_val, title_fig=False, only_half1=False, only_half2=False, cmap_used1='seismic', cmap_used2='RdBu_r'):
    """
    Plots two correlation matrices side by side with corresponding p-values.

    Args:
        df_corr1 (DataFrame): The first correlation matrix.
        df_pval1 (DataFrame): The p-values corresponding to the first correlation matrix.
        df_corr2 (DataFrame): The second correlation matrix.
        df_pval2 (DataFrame): The p-values corresponding to the second correlation matrix.
        figsize_val (tuple): The size of the figure (width, height).
        width_ratios_val (list): The width ratios of the subplots.
        title_fig (str, optional): The title of the figure. Defaults to False.
        only_half1 (bool, optional): If True, only plots the upper triangular part of the first correlation matrix. Defaults to False.
        only_half2 (bool, optional): If True, only plots the upper triangular part of the second correlation matrix. Defaults to False.
        cmap_used1 (str, optional): The colormap used for the first correlation matrix. Defaults to 'seismic'.
        cmap_used2 (str, optional): The colormap used for the second correlation matrix. Defaults to 'RdBu_r'.

    Returns:
        None
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=figsize_val, gridspec_kw={'width_ratios': width_ratios_val})

    if only_half1:
        fig1 = mask_plot_corr_pval(df_corr1, df_pval1, cmap_used1, ax1)
    
    elif not only_half1:
        fig1 = no_mask_plot_corr_pval(df_corr1, df_pval1, cmap_used1, ax1)

    fig1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
    
    if only_half2:
        fig2 = mask_plot_corr_pval(df_corr2, df_pval2, cmap_used2, ax2)

    elif not only_half2:
        fig2 = no_mask_plot_corr_pval(df_corr2, df_pval2, cmap_used2, ax2)

    fig2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')

    plt.show()

    if title_fig:
        fig.savefig(f'{DIR_OUTPUT}{title_fig}', dpi=300, format='jpg', bbox_inches='tight', pad_inches=0.1)


def simple_gcorr_plot(df_):
    """
    Plot a heatmap of the transposed dataframe `df_` with genetic correlation values.

    Parameters:
    df_ (pd.DataFrame): The input dataframe containing genetic correlation values.

    Returns:
    None
    """
    vmax_val = plot_limits(df_)
    fig, ax = plt.subplots(figsize=(18, 15))
    sns.heatmap(df_.T, 
                cbar=True, #If not False
                fmt="", #annot_kws={'weight': 'normal'}, 
                vmin=-vmax_val,
                vmax=vmax_val, 
                cmap='seismic',alpha=1.0, cbar_kws={'label': 'genetic correlation'},
                ax=ax)
    plt.xticks(rotation=45, ha='right')

    
def individual_star(df_corr_T, df_log_T, val_conversion_factor=2, val_conversion_factor_extra=1, title_fig=False, only_half=False, cmap_used='seismic'):
    """
    Plot a heatmap of genetic correlation values.

    Args:
        df_corr_T (DataFrame): The input DataFrame containing genetic correlation values.
        df_log_T (DataFrame): The input DataFrame containing annotation values.
        val_conversion_factor (int, optional): The conversion factor for adjusting the figure size based on the length of x-axis labels. Defaults to 2.
        val_conversion_factor_extra (int, optional): The conversion factor for adjusting the figure size based on the length of y-axis labels. Defaults to 1.
        title_fig (str or False, optional): The title of the figure. If False, no title is added. Defaults to False.
        only_half (bool, optional): If True, only the lower triangular part of the heatmap is shown. Defaults to False.
        cmap_used (str, optional): The colormap used for the heatmap. Defaults to 'seismic'.

    Returns:
        None
    """

    # Adjust the figure size based on the length of x-axis labels
    num_columns = get_size_val(df_corr_T, val_conversion_factor)
    num_rows = get_size_val(df_corr_T, val_conversion_factor, False)

    f = plt.figure()
    fig, ax = plt.subplots(figsize=(num_columns/val_conversion_factor_extra, num_rows))

    if only_half:
        mask = np.tril(np.ones_like(df_corr_T.corr(), dtype=bool))
        df_corr_T = df_corr_T.where(mask)

        vmax_val = plot_limits(df_corr_T)

        fig1 = sns.heatmap(df_corr_T, 
                           annot=df_log_T, 
                           cbar=True, #If not False
                           fmt="", annot_kws={'weight': 'normal'}, 
                           vmin=-vmax_val, 
                           vmax=vmax_val, 
                           cmap=cmap_used, alpha=1.0, cbar_kws={'label': 'Genetic correlation'},
                           mask=~mask  # Aplicar la máscara para mostrar solo la mitad inferior
                          )
    else:
        vmax_val = plot_limits(df_corr_T)

        fig1 = sns.heatmap(df_corr_T, 
                           annot=df_log_T, 
                           cbar=True, #If not False
                           fmt="", annot_kws={'weight': 'normal'}, 
                           vmin=-vmax_val, 
                           vmax=vmax_val,
                           cmap=cmap_used, alpha=1.0, cbar_kws={'label': 'Genetic correlation'}
                          )

    fig1.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    if title_fig != False:
        fig.savefig(f'{DIR_OUTPUT}{title_fig}', dpi=300, format='jpg', bbox_inches='tight', pad_inches=0.1)


def get_genetic_names(pheno_info_file):
    """
    Retrieves genetic names from a phenotype information file.

    Args:
        pheno_info_file (str): The path to the phenotype information file.

    Returns:
        tuple: A tuple containing two lists. The first list contains the genetic names with removed extensions,
               and the second list contains the genetic names as they are in the 'IDP_name' column of the file.
    """
    inf = pd.read_csv(pheno_info_file)
    inf_main = inf[inf['main_IDP'] == 'TRUE']

    diseases_traits_old_intermed = list(inf_main['used_file'])
    ## We need to remove the extension
    traits_all = [re.sub(r'_modified\.(tsv|csv)', '', str(file_name)) for file_name in diseases_traits_old_intermed]
    traits_all = [filename.replace('.txt', '') for filename in traits_all]
    traits_all_new = list(inf_main['IDP_name'])

    return traits_all, traits_all_new

def get_gcorr_file_names(traits_all, list_retina_homologous_red):
    """
    Retrieves the file names for genetic correlation analysis.

    Args:
        traits_all (list): A list of traits.
        list_retina_homologous_red (list): A list of retina homologous traits.

    Returns:
        DataFrame: A DataFrame containing the file names and corresponding traits.
    """
    #save_path = DIR_LDSR
    l_diseases_all=[]
    for trait in traits_all:
        #logger.info(trait)
        for file in os.listdir(DIR_LDSR):
            #logger.info(file, '\n')
            if file.startswith(trait) and ((file.endswith('.tsv'))  or (file.endswith('.txt'))or (file.endswith('.csv'))):
                #logger.info(file)
                data={
                    'pheno':  trait,
                    'file':  file
                    #,'N': df_ss['N'].iloc[0],
                    }
                l_diseases_all.append(data)

    df_diseases_all =pd.DataFrame(l_diseases_all)
    #file_name_end = '_irnt.gwas.imputed_v3.both_sexes.tsv'

    l_traits_file=[]
    for trait in list_retina_homologous_red:
        file_pheno = f'{trait}__munged.sumstats.gz'
        l_traits_file.append(file_pheno)

    #traits_files = l_traits_file + list(df_diseases_all['file'])
    #traits_names = list_retina_homologous_red + list(df_diseases_all['pheno'])
    #logger.info(len(traits_files), len(l_traits_file), len(list(df_diseases_all['file'])))
    return df_diseases_all


# filter the files names containing 2 traits
def read_ldsr(traits_files, traits_col_index):  # sourcery skip: low-code-quality
    """
    Reads LD score regression results from log files and returns dataframes for covariance, correlation, standard deviation,
    p-value, and h2 values.

    Args:
        traits_files (list): List of trait file names.
        traits_col_index (list): List of trait column indices.

    Returns:
        df_cov (DataFrame): Dataframe for covariance values.
        df_corr (DataFrame): Dataframe for correlation values.
        df_std (DataFrame): Dataframe for standard deviation values.
        df_pval (DataFrame): Dataframe for p-value values.
        df_h2 (DataFrame): Dataframe for h2 values.
        df_h2_std (DataFrame): Dataframe for h2 and std values.
    """
    df_cov = pd.DataFrame(columns=traits_col_index, index=traits_col_index)
    df_corr = pd.DataFrame(columns=traits_col_index, index=traits_col_index)
    df_std = pd.DataFrame(columns=traits_col_index, index=traits_col_index)
    df_pval = pd.DataFrame(columns=traits_col_index, index=traits_col_index)
    df_h2 = pd.DataFrame(columns=traits_col_index, index=traits_col_index)
    df_h2_std = pd.DataFrame(columns=traits_col_index, index=traits_col_index)

    for i in range(len(traits_files)):
        for j in range(len(traits_files)):
            h2 = []
            if not traits_files[i].endswith('.sumstats.gz') and not traits_files[j].endswith('.sumstats.gz'):
                file_both_name = traits_files[i] + '.sumstats.gz' + '_' + traits_files[j] + '.sumstats.gz.log'
            elif not traits_files[i].endswith('.sumstats.gz') and traits_files[j].endswith('.sumstats.gz'):
                file_both_name = traits_files[i] + '.sumstats.gz' + '_' + traits_files[j] + '.log'
            elif traits_files[i].endswith('.sumstats.gz') and not traits_files[j].endswith('.sumstats.gz'):
                file_both_name = traits_files[i] + '_' + traits_files[j] + '.sumstats.gz.log'
            else:
                file_both_name = traits_files[i] + '_' + traits_files[j] + '.log'
            dir_traitsfile = DIR_LDSR + file_both_name

            # Check if the file exists, and if not, continue to the next iteration
            if not os.path.isfile(dir_traitsfile):
                logger.warning('ERROR WITH:', dir_traitsfile)
                continue

            with open(dir_traitsfile) as fp:
                Lines = fp.readlines()
                for line in Lines:
                    split = line.split()
                    if 'gencov:' in split:
                        df_cov.iloc[i][j] = round(float(split[split.index('gencov:') + 1]), 2)
                        df_cov.iloc[j][i] = round(float(split[split.index('gencov:') + 1]), 2)
                    if 'Correlation:' in split:
                        df_corr.iloc[i][j] = round(float(split[split.index('Correlation:') + 1]), 2)
                        df_corr.iloc[j][i] = round(float(split[split.index('Correlation:') + 1]), 2)
                        df_std.iloc[i][j] = split[3]
                        df_std.iloc[j][i] = split[3]
                    if 'P:' in split:
                        df_pval.iloc[i][j] = float(split[split.index('P:') + 1])
                        df_pval.iloc[j][i] = float(split[split.index('P:') + 1])
                    if 'Total Observed scale h2:' in line:
                        df_h2.iloc[i][j] = line
                        df_h2_std.iloc[i][j] = line
                        #match = re.search(r'Total Observed scale h2:\s*([0-9.]+)', line)
                        match = re.search(r'Total Observed scale h2:\s*([0-9.]+)\s*\(([0-9.]+)\)', line)

                        if match:
                            main_value = match.group(1)
                            df_h2.iloc[i][j] = float(main_value)
                            parenthesis_value = match.group(2)
                            df_h2_std.iloc[i][j] = f"{main_value} ({parenthesis_value})"
    
    return df_cov, df_corr, df_std, df_pval, df_h2, df_h2_std


def rename_squared(df, l_phenos_old, l_phenos_new):
    """
    Renames the columns and indexes of a DataFrame.

    Args:
        df (pandas.DataFrame): The DataFrame to be modified.
        l_phenos_old (list): A list of old column/index names.
        l_phenos_new (list): A list of new column/index names.

    Returns:
        pandas.DataFrame: The modified DataFrame with renamed columns and indexes.
    """
    df.rename(columns=dict(zip(l_phenos_old, l_phenos_new)), inplace=True)
    df.rename(index=dict(zip(l_phenos_old, l_phenos_new)), inplace=True)
    return df

def convert_df_figures(df_corr_simpl, df_std_simpl):
    """
    Convert the given dataframes `df_corr_simpl` and `df_std_simpl` by performing various operations.

    Parameters:
    df_corr_simpl (DataFrame): The correlation dataframe.
    df_std_simpl (DataFrame): The standard deviation dataframe.

    Returns:
    df_corr_simpl (DataFrame): The converted correlation dataframe.
    df_corr_minus_std (DataFrame): The converted dataframe obtained by subtracting the absolute values of `df_std_simpl` from `df_corr_simpl`.
    df_std_simpl (DataFrame): The converted standard deviation dataframe.
    """
    #df_std_simpl.drop_duplicates(inplace=True)
    df_std_simpl = df_std_simpl.loc[:,~df_std_simpl.columns.duplicated()].copy()
    df_std_simpl.columns

    #df_std_simpl = df_std_simpl.astype(str)
    for col in df_std_simpl.columns:
        #logger.info(col)
        #logger.info(df_std_simpl[col].head(3))
        df_std_simpl[col] = df_std_simpl[col].str.replace("(", "", regex=True)
        df_std_simpl[col] = df_std_simpl[col].str.replace(")", "", regex=True)

    df_std_simpl=df_std_simpl.astype(float)
    df_std_simpl=df_std_simpl.round(2)
    df_std_simpl.dtypes

    #logger.info(df_corr_simpl.head(1))
    ###df_corr_simpl= -1*df_corr_simpl.copy() ######harcoded!!! just for neale

    # Define a condition to determine whether to replace NaN with 0
    replace_nan_with_0 = False  # Set this to True if you want to replace NaN with 0

    # Conditionally replace NaN values with 0
    if replace_nan_with_0:
        df_corr_simpl = df_corr_simpl.fillna(0)
        df_std_simpl = df_std_simpl.fillna(0)

    df_corr_minus_std= (np.sign(df_corr_simpl))*((abs(df_corr_simpl) - abs(df_std_simpl)))
    df_corr_minus_std = df_corr_minus_std.astype(float)

    #logger.info(df_corr_simpl.head(1))
    df = df_corr_simpl.astype(str) + '\n (' + df_std_simpl.astype(str)+ ')'

    df_corr_minus_std= df_corr_minus_std.T
    df=df.T

    df_corr_simpl = df_corr_simpl.astype(float)

    return df_corr_simpl, df_corr_minus_std, df_std_simpl

def detele_col_index(df, l_phenos_new):
    """
    Deletes columns from a DataFrame based on a list of column names and returns the modified DataFrame.

    Parameters:
    df (DataFrame): The input DataFrame.
    l_phenos_new (list): A list of column names to be deleted.

    Returns:
    DataFrame: The modified DataFrame with the specified columns removed.
    """
    df = df.drop(columns=l_phenos_new)
    df = df[df.index.isin(l_phenos_new)]
    return df

##### pascal genes and path:
def setting_genes_path(gen_path):
    """
    Sets the parameters based on the given `gen_path` value.

    Args:
        gen_path (str): The type of genes path. Can be either 'gen' or 'pathway'.

    Returns:
        tuple: A tuple containing the following parameters:
            - sep_used (str): The separator used in the data file.
            - col_name (str): The column name for the number of genes or pathways.
            - pval_min (float): The minimum p-value threshold.
            - names_used (list): The list of column names used in the data file.
    """
    if gen_path == 'gen':
        sep_used = '\t'
        col_name = 'N genes'
        pval_min = float(P_VAL_GENES)
        names_used = [gen_path, 'p']

    elif gen_path == 'pathway':
        sep_used = ' '
        col_name = 'N pathway'
        pval_min = float(P_VAL_PATHS)
        names_used = ['pathway', 'one', 'second', 'p']

    return sep_used, col_name, pval_min, names_used

def pascal_scores_reading(gen_path, dir_used, list_files):
    """
    Reads Pascal scores from multiple files and returns a DataFrame with the results.

    Args:
        gen_path (str): The column name for the gene/pathway information in the files.
        dir_used (str): The directory where the files are located.
        list_files (list): A list of file names to read.

    Returns:
        pandas.DataFrame: A DataFrame containing the processed data.

    Raises:
        None.

    """
    l_all = []
    sep_used, col_name, pval_min, names_used = setting_genes_path(gen_path)
    os.chdir(dir_used)
    for file in list_files:
        try:
            df = pd.read_csv(file, delimiter=sep_used, names=names_used)
            df['file_col'] = file
            #logger.info(file, df.head(2))
            df['-log10(p)'] = -np.log10(df['p'])
            y = df[df['-log10(p)'] >= pval_min]
            df_significant = y.sort_values('-log10(p)', ascending=False)

            if "_modified" in file:
                pheno, b = file.split("_modified.")
            else:
                pheno, b = file.split("__")

            data = {
                'name_pheno': pheno,
                col_name: len(df_significant),
                'names': df_significant[gen_path].to_list()
            }
            if gen_path == 'path':
                data['pathway_names'] = df_significant['pathway'].to_list()
            
            l_all.append(data)
            
        except Exception:
            logger.info(f'No possible to open {file}')
            continue

    return pd.DataFrame(l_all)
    

### READ OTHER VASCULAR IDPS
def get_list_files_dir(dir_used, ending_used):
    """
    Get a list of files in a directory with a specific file extension.

    Args:
        dir_used (str): The directory path.
        ending_used (str): The file extension to filter the files.

    Returns:
        list: A list of file paths that match the specified file extension.
    """
    os.chdir(dir_used)
    ending_files = ending_used
    return(glob.glob(ending_files))


def complete_genes_path(gen_path, list_retina_homologous_red, list_retina_homologous_red_new, list_names, traits_all_old, traits_all_new):
    """
    Process and complete the genes path.

    Args:
        gen_path (str): The type of genes path ('pathway' or 'gen').
        list_retina_homologous_red (list): List of retina homologous genes.
        list_retina_homologous_red_new (list): List of new retina homologous genes.
        list_names (list): List of gene names.
        traits_all_old (list): List of old traits.
        traits_all_new (list): List of new traits.

    Returns:
        tuple: A tuple containing two dataframes:
            - df_shared_elements_names: Dataframe containing shared elements and their names.
            - df_file_form_2: Dataframe with reordered and renamed columns.

    """
    if gen_path=='pathway':
        ending = "__pathway_scores.txt"
    elif gen_path=='gen':
        ending= "__gene_scores"
    else:
        logger.error(f'Not a proper ending {ending}')
    ending_used = f'*{ending}'

    ### READ OTHER VASCULAR IDPS
    list_files_found = get_list_files_dir(DIR_OTHER_GENES_PATH, ending_used)
    df_other = pascal_scores_reading(gen_path, DIR_OTHER_GENES_PATH, list_files_found)
    # Use str.replace to remove '.txt' from the 'filename' column
    df_other['name_pheno'] = df_other['name_pheno'].str.replace('.txt', '', regex=False)

    ### READ RETINA IDPS
    list_of_files = [f'{element}{ending}' for element in list_retina_homologous_red]
    df_retina = pascal_scores_reading(gen_path, DIR_GENES_RET, list_of_files)
    df_all = pd.concat([df_retina, df_other], ignore_index=True)

    df_shared_elements, df_shared_elements_names = calculate_shared_elements_counts(df_all)

    df_shared_elements_2 = df_shared_elements.T
    df_shared_elements_2 = reorder_index_columns(df_shared_elements_2, list_names)
    df_shared_elements_2 = df_shared_elements_2.replace(np.nan, 0)
    df_shared_elements_2 = df_shared_elements_2.astype(np.int64)
    df_shared_elements_2 = rename_col_index(df_shared_elements_2, traits_all_old+list_retina_homologous_red, traits_all_new+list_retina_homologous_red_new)

    df_shared_elements_names_2 = reorder_index_columns(df_shared_elements_names, list_names)
    df_shared_elements_names_2 = rename_col_index(df_shared_elements_names_2, traits_all_old+list_retina_homologous_red, traits_all_new+list_retina_homologous_red_new)

    return df_shared_elements_names_2, df_shared_elements_2


def pascal_scatter_names(gen_path, list_retina_homologous_red=False):
    """
    Generate scatter plots for Pascal dataset.

    Args:
        gen_path (str): The type of genes to consider. Can be 'pathway' or 'gen'.
        list_retina_homologous_red (bool, optional): Whether to consider retina IDPs. Defaults to False.

    Returns:
        tuple: A tuple containing the following dataframes:
            - df_significant: Dataframe of significant genes/pathways per phenotype.
            - df_count: Dataframe of the number of significant genes/pathways per phenotype.
            - df_guardar_final: Dataframe of the significant genes/pathways per phenotype per file.
            - df_save_shapes: Dataframe of the intersection lengths between different files.
            - df_save_intersections: Dataframe of the genes/pathways in the intersections between different files.
            - list_files_found: List of files found.

    Raises:
        Exception: If there is an error opening a file.
    """
    l_aux = []
    if gen_path=='pathway':
        ending = "__pathway_scores.txt"
    elif gen_path=='gen':
        ending= "__gene_scores"
    else:
        logger.error(f'Not a proper ending {ending}')
        
    ending_used = f'*{ending}'

    if list_retina_homologous_red==False:
        ### READ OTHER VASCULAR IDPS
        dir_used = DIR_OTHER_GENES_PATH
        list_files_found = get_list_files_dir(dir_used, ending_used)

    else:
        ### READ RETINA IDPS
        dir_used = DIR_GENES_RET 
        list_files_found = [f'{element}{ending}' for element in list_retina_homologous_red]
            
    sep_used, col_name, pval_min, names_used = setting_genes_path(gen_path)
    os.chdir(dir_used)
    for file in list_files_found:
        try:
            df = pd.read_csv(file, delimiter=sep_used, names=names_used)
            df['file_col'] = file
            l_aux.append(df)
            # Concat all the csvs
            df_concat = pd.concat(l_aux)
            #From p to -log10(p)
            df_concat['-log10(p)'] = -np.log10(df_concat['p'])
            y = df_concat[df_concat['-log10(p)'] >= pval_min]
            df_significant = y.sort_values('-log10(p)', ascending=False)

            df_count=df_significant[gen_path].value_counts().to_frame()
            df_count['ratio_N_pheno']=df_significant[gen_path].value_counts().to_frame()/len(list_files_found)

            ## Save the number of significant genes/path per phenotype
            df_guardar = pd.DataFrame(df_significant.groupby(by=['file_col'])[gen_path].apply(list))
            df_guardar2 = pd.DataFrame(df_significant.groupby(by=['file_col'])[gen_path].count())
            df_guardar_final=df_guardar.merge(df_guardar2, how='inner', on='file_col')
            df_save_shapes=pd.DataFrame([])
            df_save_intersections=pd.DataFrame([])

            i=0
            for file in list_files_found:
                i=i+1
                genes=df_significant[df_significant['file_col']==file][gen_path]#.to_list()
                genes=genes.to_list()
                l_aux2 = []
                l_aux3 = []

                for j in range(len(list_files_found)):#-i):#-file: #Error
                    other_file=list_files_found[j]
                    df_intersection=df_significant[(df_significant['file_col']==other_file)&(df_significant[gen_path].isin(genes))]
                    # To save the intersection len
                    save_shapes=df_intersection.shape[0]
                    l_aux3.append(save_shapes)

                    # To save the names of the genes in the intersection
                    l_aux2.append(df_intersection[gen_path].to_list())# {'index1': value1, 'index2':value2,...}, ignore_index=True)#.to_list())
                    #logger.info('antes', df_intersection[gen_path].to_list(), 'despues')
                    #logger.info(df_intersection[gen_path].values())

                # To save the intersection len
                df = pd.DataFrame({file:l_aux3})
                df_save_shapes = pd.concat([df_save_shapes, df], axis=1)

                # To save the names of the genes in the intersection
                df2 = pd.DataFrame({file:l_aux2})
                df_save_intersections = pd.concat([df_save_intersections, df2], axis=1)
                #logger.info(df_save_intersections[1])

            # To save the intersection len  
            df_save_shapes = df_save_shapes.set_axis(df_save_shapes.columns, axis='index')

            df_save_intersections = df_save_intersections.set_axis(df_save_intersections.columns, axis='index')
        
        except Exception:
            logger.warning('No possible to open:', file)
            continue
        
    return df_significant, df_count, df_guardar_final, df_save_shapes, df_save_intersections, list_files_found


def combine_genes_path_names_idps_ret(df1, df2):
    """
    Combines two DataFrames vertically and returns the combined DataFrame.

    Args:
        df1 (pandas.DataFrame): The first DataFrame.
        df2 (pandas.DataFrame): The second DataFrame.

    Returns:
        pandas.DataFrame: The combined DataFrame.

    Raises:
        AssertionError: If the shapes of the DataFrames don't match after concatenation.
    """
    # Check sizes of original DataFrames
    logger.info(f"Size of df1: {df1.shape}. Size of df2: {df2.shape}")

    # Check size of combined DataFrame
    combined_df = pd.concat([df1, df2], ignore_index=True)
    logger.info(f"Size of combined_df is {combined_df.shape}")

    assert combined_df.shape == (df1.shape[0] + df2.shape[0], df1.shape[1]), "Error: The shapes of the DataFrames don't match after concatenation."
    
    return combined_df
    
def reorder_yaxis_gen_path_names(df_reduc, gen_path):
    """
    Reorders the y-axis of a DataFrame based on the count of points for each gen.

    Args:
        df_reduc (pandas.DataFrame): The DataFrame to be reordered.
        gen_path (str): The column name representing the gen path.

    Returns:
        pandas.DataFrame: The reordered DataFrame.
    """
    ### RE ORDER THE Y AXIS
    # Calculate the count of points for each gen
    gen_counts = df_reduc[gen_path].value_counts()

    # Sort the gen values based on the count
    sorted_gens = gen_counts.sort_values(ascending=False).index

    # Create a new column to specify the order of the gens
    df_reduc['gen_order'] = df_reduc[gen_path].map(dict(zip(sorted_gens, range(len(sorted_gens)))))

    ### not completely sorted, but the best solution I found
    df_reduc.sort_values(by=['order', 'gen_order'], ascending=[True, True], inplace=True)

    return df_reduc

def df_gen_path_names_processing(all_names, all_names_new, df_significant, N_head, gen_path):
    """
    Process the DataFrame `df_significant` by filtering and reordering the data based on the given parameters.

    Args:
        all_names (list): List of all original names.
        all_names_new (list): List of corresponding new names.
        df_significant (pandas.DataFrame): The input DataFrame to be processed.
        N_head (int): Number of top genes/pathways to include in the output DataFrame.
        gen_path (str): Type of path (either 'gen' or 'pathway').

    Returns:
        pandas.DataFrame: The processed DataFrame with filtered and reordered data.
    """
    
    if gen_path == 'gen':
        csv_regex = '_modified.csv__gene_scores'
        tsv_regex = '_modified.tsv__gene_scores'
        txt_regex = '.txt__gene_scores'
        ret_regex = '__gene_scores'
        #set_no_asso = {'1365'}

    elif gen_path == 'pathway':
        csv_regex = '_modified.csv__pathway_scores.txt'
        tsv_regex = '_modified.tsv__pathway_scores.txt'
        txt_regex = '.txt__pathway_scores.txt'
        ret_regex = '__pathway_scores.txt'
        #using brain vessels I and vol
        #set_no_asso = {'medianDiameter_artery', '1365', '0220', 'medianDiameter_vein', 'AAdis_model1_bolt_P_BOLT_LMM_INF', 'eq_CRAE', '0203', '1351', 'DAdis_model1_bolt_P_BOLT_LMM_INF'}
        #set_no_asso = {'medianDiameter_artery', 'medianDiameter_vein', 'AAdis_model1_bolt_P_BOLT_LMM_INF', 'eq_CRAE', '0203', '1351', 'DAdis_model1_bolt_P_BOLT_LMM_INF', 'deep_adj_vol'} #, 'pvent_adj_vol', 'deep_adj_vol'}

    
    # Create mapping dictionaries
    mapping_labels = dict(zip(all_names, all_names_new))
    mapping_order = dict(zip(all_names, list(range(len(all_names)))))

    # Copy and filter DataFrame
    df_2 = df_significant.copy()
    df_2['file_col'] = df_2['file_col'].str.replace(csv_regex, '',  regex=True)
    df_2['file_col'] = df_2['file_col'].str.replace(tsv_regex, '',  regex=True)
    df_2['file_col'] = df_2['file_col'].str.replace(txt_regex, '',  regex=True)
    df_2['file_col'] = df_2['file_col'].str.replace(ret_regex, '',  regex=True)
    #### SOME DO NOT HAVE VALUES, BECAUSE THEY HAVE NOT SIGN GENES/PATH
    df_i = df_2[df_2['file_col'].isin(all_names)]

    # Assert the condition
    unique_names_df1 = set(df_i['file_col'])
    unique_names_df2 = set(all_names)
    #logger.info(len(unique_names_df1), len(all_names))
    # Find names present in df1 but not in df2
    names_only_in_df1 = unique_names_df1 - unique_names_df2
    names_only_in_df2 = unique_names_df2 - unique_names_df1 
    #logger.info('FIRST', len(names_only_in_df1),names_only_in_df1)
    #logger.info('SECOND', len(names_only_in_df2),names_only_in_df2)
    # if not names_only_in_df1.issubset(set_no_asso) or not names_only_in_df2.issubset(set_no_asso):
    #     raise ValueError('ERROR! The phenotypes without associations are NOT the expected ones!')
    # else:
    #     logger.info('Good. The phenotypes without associations are the expected ones!')
        
    list_genes = list(df_i[gen_path].value_counts().head(N_head).index)

    df_reduc = df_i[df_i[gen_path].isin(list_genes)].copy()

    df_reduc['file_col'] = df_reduc['file_col'].map(mapping_labels)
    df_reduc['order']= df_reduc['file_col'].map(mapping_order)

    df_reduc.sort_values(by='order', ascending=True, inplace=True)

    df_reduc = reorder_yaxis_gen_path_names(df_reduc, gen_path)

    return df_reduc

def plot_scatter_pascal_names(df_reduc, gen_path, fig_name):
    """
    Plot scatter plot with variable-sized markers based on -log10(p) values.

    Args:
        df_reduc (DataFrame): The reduced DataFrame containing the data.
        gen_path (str): The column name for the y-axis values.
        fig_name (str): The name of the output figure file.

    Returns:
        None
    """
    cte = 0.5
    x_size = len(df_reduc['file_col'].unique())
    y_size = len(df_reduc[gen_path].unique())

    figsize_val = (cte * x_size, cte * y_size)

    fig, ax = plt.subplots(1, 1, figsize=figsize_val)

    fig2 = ax.scatter(df_reduc['file_col'], df_reduc[gen_path], s=1.25 * df_reduc['-log10(p)'], c='grey')
    plt.grid(color='gray', linestyle='-', linewidth=0.1)

    cmap_label = 'Pleiotropy (N trait pairs)'
    plt.xticks(rotation=45, ha='right')
    kw = dict(prop='sizes', num=4, color='grey', alpha=0.6)

    ax.legend(*fig2.legend_elements(**kw), loc='upper right', title='Mean -log10(p)', bbox_to_anchor=(1, 1.25), title_fontsize=12)

    fig.savefig(f'{DIR_OUTPUT}{fig_name}', dpi=300, format='jpg', bbox_inches='tight', pad_inches=0.1)


def read_and_process_csv(file_path):
    df = pd.read_csv(file_path)
    """
    Reads a CSV file from the given file path and processes it.

    Parameters:
    file_path (str): The path to the CSV file.

    Returns:
    pandas.DataFrame: The processed DataFrame.
    """
    df.index = df['Unnamed: 0']
    df.drop(columns=['Unnamed: 0'], inplace=True)
    df.index.name = None
    return df

def read_and_process_log10p_csv(file_path):
    """
    Reads a CSV file from the given file path and processes it.

    Args:
        file_path (str): The path to the CSV file.

    Returns:
        pandas.DataFrame: The processed DataFrame.
    """
    df = pd.read_csv(file_path)
    df.replace(np.nan, '', inplace=True)
    df.index = df['Unnamed: 0']
    df.drop(columns=['Unnamed: 0'], inplace=True)
    df.index.name = None
    return df

def squared_detele_col_index(df, l_phenos_new):
    """
    Returns a DataFrame with only the columns specified in l_phenos_new.

    Parameters:
    df (DataFrame): The input DataFrame.
    l_phenos_new (list): A list of column names to keep in the DataFrame.

    Returns:
    DataFrame: A new DataFrame with only the specified columns.
    """
    df = df[l_phenos_new]
    return df[df.index.isin(l_phenos_new)]
