import matplotlib.pyplot as plt

def configure_global_params():
    plt.rcParams['font.size'] = 16
    plt.rcParams['figure.facecolor'] = 'white'

FIG_PHENO_GENO_SIZE = (12, 2.5)
FIG_PHENO_GENO_RATIOS = [0.6, 1]

def pheno_img_values(IDP_VASCULAR_USED, combine_IMT_heart_brain, main_final_figs):
    """
    Calculate and return the values for various plotting settings based on the input parameters.

    Parameters:
    - IDP_VASCULAR_USED (str): The type of vascular used. Should be either 'homologous_all' or 'all'.
    - combine_IMT_heart_brain (bool): Whether to combine IMT for heart and brain.
    - main_final_figs (bool): Whether the figures are for main final plots.

    Returns:
    - figsize_val_square (tuple): The figure size for square plots.
    - cbar_fraction_val_square (float): The fraction of colorbar size for square plots.
    - cte_square (float): A constant value for square plots.
    - figsize_val_both (tuple): The figure size for combined plots.
    - width_ratios_val_both (list): The width ratios for combined plots.
    """
    if IDP_VASCULAR_USED not in ['homologous_all', 'all']:
        return
    if (combine_IMT_heart_brain == False) and (main_final_figs == False): ## Sup figs 
        figsize_val_square = (2*54,54)
        cbar_fraction_val_square = 0.045
        cte_square = 4.5
        figsize_val_both = (30, 20)
        width_ratios_val_both = [0.8, 0.2]

    if combine_IMT_heart_brain == True:
        if main_final_figs == False:
            figsize_val_square = (2*24,24)
            cbar_fraction_val_square = 0.04
            cte_square = 3
            figsize_val_both = (32.4, 18)
            width_ratios_val_both = [0.8, 0.25]

        elif main_final_figs == True: #main
            figsize_val_square = (16,16)
            cbar_fraction_val_square = 0.04
            cte_square = 3
            figsize_val_both = FIG_PHENO_GENO_SIZE #(12, 4)
            width_ratios_val_both = FIG_PHENO_GENO_RATIOS #if cols 5, 11 #[1, 1] #if cols 10, 11

    return figsize_val_square, cbar_fraction_val_square, cte_square, figsize_val_both, width_ratios_val_both
