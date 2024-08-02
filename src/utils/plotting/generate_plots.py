import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr#, kendalltau, spearmanr
from matplotlib.colors import LogNorm
import seaborn as sns
import logging

from utils.plotting.plotting_settings import configure_global_params
from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS, DIR_OUTPUT, DATE_USED, ALPHA_1, ALPHA_2

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)

configure_global_params()

sys.path.append(DIR_UTILS)

#def kendall_pval(x,y) d:
#    return kendalltau(x,y)[1]
#def spearmanr_pval(x,y):
#    return spearmanr(x,y)[1]

def pearsonr_pval(x,y):
    return pearsonr(x,y)[1]

def plot_cov_histograms(df_cov, file_name):
    """
    Plots histograms for each column in the given DataFrame.

    Args:
        df_cov (pandas.DataFrame): The DataFrame containing the data.
        file_name (str): The name of the file to save the plot.

    Returns:
        None
    """
    num_columns = df_cov.shape[1]

    # Compute the num of rows and cols needed for the subplots
    num_rows = (num_columns - 1) // 4 + 1  # 4 subplots per row

    fig, axes = plt.subplots(num_rows, 4, figsize=(20, 5 * num_rows), sharey=True)
    axes = axes.flatten()

    for i, column in enumerate(df_cov.columns):
        if i < num_columns:
            try:
                axes[i].hist(df_cov[column], bins=10)
            except RuntimeWarning as e:
                print(f"Warning: {e}")
                print(f"Skipping histogram for column: {column}")
                continue
            axes[i].set_title(column)
            # Calculate the number of non-null values in the column
            not_null_count = df_cov[column].count()
            axes[i].legend([f"Not Null: {not_null_count}"])
    plt.tight_layout() # Adjust space bewtween the subplots
    save_path = os.path.join(DIR_OUTPUT, file_name)
    print(f'Saving cov histograms... {save_path}')
    plt.savefig(save_path)
    plt.show()

def filtered_corr_plots(df, dict_subset, figsize_1, cbar_1, SAVE_FIGURES, name_used, only_half):
    """
    Generate filtered correlation plots for a given DataFrame.

    Parameters:
    - df (pandas.DataFrame): The input DataFrame.
    - dict_subset (dict): A dictionary containing the column names to be included in the correlation analysis.
    - figsize_1 (tuple): A tuple specifying the figure size of the correlation plot.
    - cbar_1 (float): The fraction of the colorbar to be displayed.
    - SAVE_FIGURES (bool): Whether to save the generated figures.
    - name_used (str): The name to be used when saving the figures.
    - only_half (bool): Whether to display only the lower triangular part of the correlation matrix.

    Returns:
    None
    """
    list_subset = list(dict_subset.values())
    df_filtered = df[list_subset]
    plot_custom_corr_simplified(df_filtered.corr(), "", figsize_1, cmap_used='seismic', ax=None, cbar_fraction=cbar_1, SAVE_FIGURES=SAVE_FIGURES, name_used=name_used, only_half=only_half)

def plot_custom_corr_simplified(df, title_name, figsize_val, cmap_used=None, ax=None, cbar_fraction=0.05, SAVE_FIGURES=False, name_used=False,only_half=False):
    """
    Plot a custom correlation matrix.

    Parameters:
    - df: pandas DataFrame
        The input DataFrame containing the correlation matrix.
    - title_name: str
        The title of the plot.
    - figsize_val: tuple
        The size of the figure (width, height) in inches.
    - cmap_used: str, optional
        The colormap to be used for the plot. Default is None.
    - ax: matplotlib Axes, optional
        The Axes object to plot on. Default is None.
    - cbar_fraction: float, optional
        The fraction of the figure width that the colorbar occupies. Default is 0.05.
    - SAVE_FIGURES: bool, optional
        Whether to save the plot as an image file. Default is False.
    - name_used: bool, optional
        Whether to include a name in the saved image file. Default is False.
    - only_half: bool, optional
        Whether to plot only the lower triangular part of the correlation matrix. Default is False.

    Returns:
    - None
    """
    if cmap_used is None:
        cmap_used = 'seismic' if np.any(df.values < 0) else 'Blues'

    # Create a sub-DataFrame with non-null values
    sub_df = df.dropna(axis=0, how='all').dropna(axis=1, how='all')

    if only_half:
        mask = np.tril(np.ones_like(sub_df.corr(), dtype=bool))
        sub_df = sub_df.where(mask)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize_val, constrained_layout=True)
    else:
        fig = ax.get_figure()

    if cmap_used == 'Blues':
        im = ax.imshow(sub_df, cmap=cmap_used, interpolation='none', norm=LogNorm())
        cbar = fig.colorbar(im, ax=ax, norm=LogNorm(), fraction=cbar_fraction)
    elif cmap_used == 'Greys':
        im = ax.imshow(sub_df, cmap=cmap_used, interpolation='none', norm=LogNorm())
        cbar = fig.colorbar(im, ax=ax, norm=LogNorm(), fraction=cbar_fraction)
        
    else:  #if title_name == 'Correlation Plot':
        im = ax.imshow(sub_df, cmap=cmap_used, interpolation='none', vmin=-1, vmax=1)  # Set colorbar limits to -1 and 1 
        cbar = fig.colorbar(im, ax=ax, fraction=cbar_fraction) 
    # Rotate x-axis labels to 45 degrees
    ax.set_xticks(range(len(sub_df.columns)))
    ax.set_xticklabels(sub_df.columns, rotation=45, ha='right')
    ax.set_yticks(range(len(sub_df.index)))
    ax.set_yticklabels(sub_df.index)

    #ax.set_title(title_name)
    # Adjust the position of the colorbar
    #cbar.ax.yaxis.set_label_position('left')
    #cbar.set_label('Colorbar Label')
    if SAVE_FIGURES:
        save_path = f'{DIR_OUTPUT}{DATE_USED}_{name_used}_correlation.jpg'
        plt.savefig(save_path, bbox_inches='tight')
    
    if ax is None:
        plt.show()

def pval_asterisk(df, df_pval, multiple_testing, soft_multiple_testing, list_values, list_retina_new, log_pval=False, N_shape=False):
    """
    Applies asterisks to p-values based on significance thresholds.

    Args:
        df (DataFrame): The input DataFrame.
        df_pval (DataFrame): The p-value DataFrame.
        multiple_testing (str): The type of multiple testing correction.
        soft_multiple_testing (str): Whether to apply soft multiple testing correction.
        list_values (list): The list of values.
        list_retina_new (list): The list of new retina values.
        log_pval (bool, optional): Whether to use logarithmic p-values. Defaults to False.
        N_shape (bool, optional): The shape of N. Defaults to False.

    Returns:
        DataFrame: The DataFrame with asterisks applied to p-values.
    """
    if N_shape == False:
        N_shape = len(df.columns) ** 2 if multiple_testing == 'Complete' else len(list_values)*(len(list_values)/2+len(list_retina_new)) if multiple_testing == 'Subset' else 1

        if multiple_testing == 'Subset' and soft_multiple_testing == 'True':
            N_shape = 2*len(list_retina_new) + len(list_values)
    N_shape = int(N_shape)
    print('N_shape', N_shape)
    if log_pval:
        Bonf_thresh = -np.log10(ALPHA_1 / N_shape)
        Bonf_thresh2 = -np.log10(ALPHA_2 / N_shape)
    else:
        Bonf_thresh = (ALPHA_1 / N_shape)
        Bonf_thresh2 = (ALPHA_2 / N_shape)

    df_pval = df_pval
    p_copy = df_pval.copy()
    p_copy2 = df_pval.copy()

    p_copy = (p_copy < Bonf_thresh).replace({True:'*', False:''})
    p_copy2 = (p_copy2 < Bonf_thresh2).replace({True:'*', False:''})
    return p_copy + p_copy2


def plot_corr_pval(df_corr, df_pval, val_conversion_factor=2, val_conversion_factor_extra=1, title_fig=False, only_half=False, cmap_used='seismic'):
    """
    Plot correlation and p-value matrices.

    Args:
        df_corr (DataFrame): The correlation matrix.
        df_pval (DataFrame): The p-value matrix.
        val_conversion_factor (int, optional): The conversion factor for adjusting the figure size based on the length of x-axis labels. Defaults to 2.
        val_conversion_factor_extra (int, optional): The extra conversion factor for adjusting the figure size. Defaults to 1.
        title_fig (str or False, optional): The title of the figure. If False, no title is added. Defaults to False.
        only_half (bool, optional): Whether to plot only the upper half of the matrices. Defaults to False.
        cmap_used (str, optional): The colormap to use for plotting. Defaults to 'seismic'.

    Returns:
        Figure: The generated figure object.
    """
    # Adjust the figure size based on the length of x-axis labels
    num_columns = get_size_val(df_corr, val_conversion_factor)
    num_rows = get_size_val(df_pval, val_conversion_factor, False)

    fig, ax = plt.subplots(figsize=(num_columns/val_conversion_factor_extra, num_rows))

    if only_half:
        fig1 = mask_plot_corr_pval(df_corr, df_pval, cmap_used, ax)
    else:    
        fig1 = no_mask_plot_corr_pval(df_corr, df_pval, cmap_used, ax)
    
    fig1.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    if title_fig != False:
        fig.savefig(f'{DIR_OUTPUT}{title_fig}.jpg', dpi=300, format='jpg', bbox_inches='tight', pad_inches=0.1)

    return fig1


def get_size_val(df, conversion_factor, col=True):
    """
    Calculate the size value of a DataFrame based on a conversion factor.

    Parameters:
    - df (pandas.DataFrame): The DataFrame for which the size value needs to be calculated.
    - conversion_factor (float): The conversion factor to be used in the calculation.
    - col (bool, optional): Determines whether to calculate the size value based on the number of columns (True) or rows (False). Default is True.

    Returns:
    - float: The size value of the DataFrame divided by the conversion factor.
    """
    num = 1 if col == False else 0
    return df.shape[num] / conversion_factor

def mask_plot_corr_pval(df_corr, df_pval, cmap_used, ax, label_cbar=None):
    """
    Plot a correlation matrix with p-values masked in the lower triangle.

    Parameters:
    - df_corr (DataFrame): The correlation matrix.
    - df_pval (DataFrame): The p-value matrix.
    - cmap_used (str or Colormap): The colormap to use for the heatmap.
    - ax (Axes): The matplotlib Axes object to draw the heatmap on.
    - label_cbar (str, optional): The label for the colorbar. Defaults to None.

    Returns:
    - ax (Axes): The matplotlib Axes object with the heatmap.

    """
    mask = np.tril(np.ones_like(df_corr.corr(), dtype=bool))
    df_corr = df_corr.where(mask)
    return sns.heatmap(df_corr,  
                annot=df_pval, 
                cbar=True,
                fmt="", 
                annot_kws={'weight': 'bold'}, 
                vmin=-abs(df_corr).max().max(), 
                vmax=abs(df_corr).max().max(), 
                cmap=cmap_used, 
                alpha=1.0, 
                cbar_kws={'label': label_cbar},
                mask=~mask,  # Apply the mask to show only the lower triangle
                ax=ax)

def doble_plot_corr_pval(df_corr1, df_pval1, df_corr2, df_pval2, figsize_val, width_ratios_val, title_fig=False, only_half1=False, only_half2=False, cmap_used1='seismic', cmap_used2='RdBu_r', cbar_BP=False, cmax_used=False):
    """
    Generate a double plot of correlation matrices with p-values.

    Args:
        df_corr1 (pandas.DataFrame): The first correlation matrix.
        df_pval1 (pandas.DataFrame): The p-values corresponding to the first correlation matrix.
        df_corr2 (pandas.DataFrame): The second correlation matrix.
        df_pval2 (pandas.DataFrame): The p-values corresponding to the second correlation matrix.
        figsize_val (tuple): The size of the figure (width, height).
        width_ratios_val (tuple): The width ratios of the subplots.
        title_fig (str, optional): The title of the saved figure. Defaults to False.
        only_half1 (bool, optional): Whether to plot only the upper triangular part of the first correlation matrix. Defaults to False.
        only_half2 (bool, optional): Whether to plot only the upper triangular part of the second correlation matrix. Defaults to False.
        cmap_used1 (str, optional): The colormap used for the first correlation matrix. Defaults to 'seismic'.
        cmap_used2 (str, optional): The colormap used for the second correlation matrix. Defaults to 'RdBu_r'.
        cmax_used (bool, optional): Whether to use a common maximum value for the color scale of both correlation matrices. Defaults to False.

    Returns:
        None
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=figsize_val, gridspec_kw={'width_ratios': width_ratios_val})

    if (df_pval1 is None) and (df_pval2 is None):
        if only_half1:
            fig1 = mask_plot_corr(df_corr1, cmap_used1, ax1)
        else:
            fig1 = no_mask_plot_corr(df_corr1, cmap_used1, ax1)
        if only_half2:
            fig2 = mask_plot_corr(df_corr2, cmap_used2, ax2)
        else:
            fig2 = no_mask_plot_corr(df_corr2, cmap_used2, ax2, other_vmax=cmax_used)
    else:
        if only_half1:
            fig1 = mask_plot_corr_pval(df_corr1, df_pval1, cmap_used1, ax1)
        else:
            fig1 = no_mask_plot_corr_pval(df_corr1, df_pval1, cmap_used1, ax1, cbar_BP)
        if only_half2:
            fig2 = mask_plot_corr_pval(df_corr2, df_pval2, cmap_used2, ax2)
        else:
            fig2 = no_mask_plot_corr_pval(df_corr2, df_pval2, cmap_used2, ax2, cbar_BP)

    fig1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
    fig2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')
    plt.show()

    if title_fig:
        fig.savefig(f'{DIR_OUTPUT}{title_fig}', dpi=300, format='jpg', bbox_inches='tight', pad_inches=0.1)

def no_mask_plot_corr_pval(df_corr, df_pval, cmap_used, ax, cbar_BP=False,label_cbar=None):
    """
    Plot a heatmap of correlation values with p-values as annotations.

    Parameters:
    - df_corr (DataFrame): The correlation values DataFrame.
    - df_pval (DataFrame): The p-values DataFrame.
    - cmap_used (str or Colormap): The colormap to use for the heatmap.
    - ax (Axes): The matplotlib Axes object to plot on.
    - label_cbar (str, optional): The label for the colorbar. Defaults to None.

    Returns:
    - ax (Axes): The matplotlib Axes object with the heatmap plot.

    """
    if cbar_BP=='Hypertense':
        v_value = round(0.06352537913355354,2)
    else:
        v_value = round(abs(df_corr).max().max(),2)

    return sns.heatmap(df_corr,  
                    annot=df_pval, #.T, 
                    cbar=True, #If not False
                    fmt="", 
                    annot_kws={'weight': 'bold'}, 
                    vmin=-v_value,  #0.06352537913355354
                    vmax=v_value, 
                    cmap=cmap_used,
                    alpha=1.0, 
                    cbar_kws={'label': label_cbar},
                    ax=ax)

def mask_plot_corr(df_corr, cmap_used, ax, label_cbar=None):
    mask = np.tril(np.ones_like(df_corr.corr(), dtype=bool))
    df_corr = df_corr.where(mask)
    return sns.heatmap(df_corr, 
                cbar=True,
                fmt="", 
                annot_kws={'weight': 'bold'}, 
                vmin=0, 
                vmax=abs(df_corr).max().max(), 
                cmap=cmap_used, 
                alpha=1.0, 
                cbar_kws={'label': label_cbar},
                mask=~mask,  # Apply the mask to only show a half of the matrix
                ax=ax)

def no_mask_plot_corr(df_corr, cmap_used, ax, label_cbar=None, other_vmax=False):
    """
    Plot a correlation heatmap without a mask.

    Parameters:
    - df_corr: pandas DataFrame
        The correlation matrix to be plotted.
    - cmap_used: str or colormap object
        The colormap to be used for the heatmap.
    - ax: matplotlib Axes object
        The Axes object where the heatmap will be plotted.
    - label_cbar: str, optional
        The label for the colorbar. Default is None.
    - other_vmax: bool or float, optional
        If False, the maximum value for the colorbar will be determined automatically based on the correlation matrix.
        If a float is provided, it will be used as the maximum value for the colorbar. Default is False.

    Returns:
    - seaborn.heatmap object
        The heatmap plot.

    """
    vmax_val = abs(df_corr).max().max() if other_vmax == False else other_vmax
    return sns.heatmap(df_corr,  
                       cbar=True,  # If not False
                       fmt="", 
                       annot_kws={'weight': 'bold'}, 
                       vmin=0, 
                       vmax=vmax_val, 
                       cmap=cmap_used,
                       alpha=1.0, 
                       cbar_kws={'label': label_cbar},
                       ax=ax)

def plot_multiple_fig_shared(df_corr_matrix, df_non_null_values_count, figsize_val=(10, 6), cbar_fraction_val=0.05):
    """
    Plot multiple figures with shared y-axis.

    Parameters:
    - df_corr_matrix (DataFrame): The correlation matrix DataFrame.
    - df_non_null_values_count (DataFrame): The DataFrame containing the count of non-null values.
    - figsize_val (tuple, optional): The size of the figure. Defaults to (10, 6).
    - cbar_fraction_val (float, optional): The fraction of the colorbar. Defaults to 0.05.
    """
    cmap_corr = 'seismic' if np.any(df_corr_matrix.values < 0) else 'Blues'
    if np.any(df_non_null_values_count.values < 0):
        cmap_values = 'seismic'
    else:
        cmap_values = 'Blues'

    # Create subplots with 1 row and 2 columns
    fig, axes = plt.subplots(1, 2, figsize=figsize_val)

    # Plot correlation matrix
    plot_custom_corr_simplified(df_corr_matrix, "Correlation Plot", figsize_val, cmap_used=cmap_corr, ax=axes[0], cbar_fraction=cbar_fraction_val)

    # Plot non-null values count matrix
    plot_custom_corr_simplified(df_non_null_values_count, "Values Plot", figsize_val, cmap_used=cmap_values, ax=axes[1], cbar_fraction=cbar_fraction_val)
    # Hide y-axis of the second plot
    axes[1].set_yticks([])
    plt.show()

######
# Filter by specific rows (index) and columns
def filter_dataframe(df, rows, columns):
    """
    Filters a dataframe based on the specified rows and columns.

    Args:
        df (pandas.DataFrame): The input dataframe.
        rows (list or slice): The rows to include in the filtered dataframe.
        columns (list or slice): The columns to include in the filtered dataframe.

    Returns:
        pandas.DataFrame: The filtered dataframe.

    """
    return df.loc[rows, columns]

def get_files_names(IDPs_retina_used, IDPs_vascular_used, z_score_applied, combine_IMT_heart_brain, covariants, DATE_USED):
    """
    Generates file names based on input parameters.

    Args:
        IDPs_retina_used (str): The IDPs used for the retina.
        IDPs_vascular_used (str): The IDPs used for the vascular system.
        z_score_applied (bool): Indicates whether z-score is applied.
        combine_IMT_heart_brain (bool): Indicates whether IMT for heart and brain are combined.
        covariants (str): The covariants used.
        DATE_USED (str): The date used.

    Returns:
        tuple: A tuple containing two file names. The first file name represents the count of values and subjects, 
               while the second file name represents the z-scored data.

    """
    common_part = f"{IDPs_retina_used}_{IDPs_vascular_used}_zcored_{z_score_applied}_combine_{combine_IMT_heart_brain}_cov_used_{covariants}.csv"
    name_file_values_count = f"{DATE_USED}_N_subjects_{common_part}"
    name_file_zcored = f"{DATE_USED}_{common_part}"
    return name_file_values_count, name_file_zcored


def get_names_figs(DATE_USED, SQUARE_FIG, DOBLE_FIG, z_score_applied, combine_IMT_heart_brain, covariants, multiple_testing, SAVE_FIGURES, type_BP, FILTER_OUTLIERS, TYPE_OF_FIGS):
    """
    Generate titles for square and double figures based on the provided parameters.

    Args:
        DATE_USED (str): The date used for the titles.
        SQUARE_FIG (str): The square figure name.
        DOBLE_FIG (str): The double figure name.
        z_score_applied (bool): Indicates whether z-score is applied.
        combine_IMT_heart_brain (bool): Indicates whether IMT heart and brain are combined.
        covariants (str): The covariants used.
        multiple_testing (str): The type of multiple testing used.
        SAVE_FIGURES (bool): Indicates whether to save the figures.
        type_BP (str): The type of BP used.
        TYPE_OF_FIGS (str): The type of figures.

    Returns:
        tuple: A tuple containing the titles for square figure, double figure, and double figure with sample size.
    """
    if SAVE_FIGURES:
        if type_BP in [False, "False"]:
            common_part = f"_phe_zcored_{z_score_applied}_combine_{combine_IMT_heart_brain}_cov_used_{covariants}_multest_{multiple_testing}_{TYPE_OF_FIGS}.jpg"
        else:
            common_part = f"_phe_zcored_{z_score_applied}_combine_{combine_IMT_heart_brain}_cov_used_{type_BP}_multest_{multiple_testing}_{TYPE_OF_FIGS}.jpg"
        if FILTER_OUTLIERS not in {"False", False}:
            common_part = f"_outliers_removed_{FILTER_OUTLIERS}_{common_part}"
        title_square = f"{DATE_USED}_{SQUARE_FIG}_{common_part}"
        title_doble = f"{DATE_USED}_{DOBLE_FIG}_{common_part}"
        title_doble_N = f"{DATE_USED}_{DOBLE_FIG}_sample_size_{TYPE_OF_FIGS}.jpg"

    else:
        title_square, title_doble, title_doble_N = False, False, False

    return title_square, title_doble, title_doble_N