import sys
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.colors import LogNorm
from sklearn.metrics import mean_squared_error, r2_score
import scipy.stats as stats
import seaborn as sns
import logging

from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS, DIR_OUTPUT, DATE_USED
from utils.plotting.plotting_settings import configure_global_params

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)

configure_global_params()

sys.path.append(DIR_UTILS)


def multiple_scatterplots(df, list_param_a, list_param_b, title):
    """
    Generate multiple scatter plots based on the given parameters.

    Parameters:
    - df (pandas.DataFrame): The DataFrame containing the data.
    - list_param_a (list): List of parameter names for the x-axis of each scatter plot.
    - list_param_b (list): List of parameter names for the y-axis of each scatter plot.
    - title (str): The title for the scatter plots.

    Returns:
    None
    """
    # Calculate the number of rows and columns for the subplot grid
    num_plots = len(list_param_a)
    num_cols = 3  # Number of columns in the grid
    num_rows = -(-num_plots // num_cols)  # Calculate the number of rows needed (ceiling division)

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(15, 5*num_rows))  # Create subplot grid

    for i, (param_a, param_b) in enumerate(zip(list_param_a, list_param_b)):
        row = i // num_cols  # Calculate the current row index
        col = i % num_cols   # Calculate the current column index

        ax = axs[row, col] if num_rows > 1 else axs[col]  # Select current subplot
        ax.scatter(df[param_a], df[param_b], alpha=0.5)   # Scatter plot

        ax.set_xlabel(param_a)
        ax.set_ylabel(param_b)
        #ax.set_title(f'cov {str(covariants)} {param_a} vs {param_b} - N={len(df)}', fontsize=10)
        ax.set_title(f'cov {title}, N={len(df[[param_a, param_b]].dropna())}')

        ax.grid(True)

    plt.tight_layout()
    plt.show()


def plot_with_trendlines(ax, df, column_x, column_y, title):
    """
    Plot a scatter plot with trendlines on the given axes.

    Parameters:
    - ax (matplotlib.axes.Axes): The axes on which to plot the data.
    - df (pandas.DataFrame): The DataFrame containing the data.
    - column_x (str): The name of the column to use as the x-axis.
    - column_y (str): The name of the column to use as the y-axis.
    - title (str): The title of the plot.

    Returns:
    None
    """
    # Scatter plot
    ax.scatter(df[column_x], df[column_y], alpha=0.5, label='Data Points')

    # Linear trendline
    linear_coefficients = np.polyfit(df[column_x], df[column_y], 1)
    linear_trendline = np.poly1d(linear_coefficients)
    ax.plot(df[column_x], linear_trendline(df[column_x]), label=f'Linear Trendline: {linear_trendline}', color='red')

    # Set labels and title
    ax.set_xlabel(column_x)
    ax.set_ylabel(column_y)
    ax.set_title(f'cov {title} - N={len(df)}')
    ax.legend()  # Show legend
    ax.grid(True)  # Add grid

def simple_scatter_plot(df, column_x, column_y, title):
    """_summary_

    Args:
        df (_type_): _description_
        column_x (_type_): _description_
        column_y (_type_): _description_
        title (_type_): _description_
    """
    # Create a scatter plot
    plt.figure(figsize=(8, 6))  # Define the size of the plot (optional)
    plt.scatter(df[column_x], df[column_y], alpha=0.5)  # Create scatter plot with transparency

    # Set labels and title
    plt.xlabel(column_x)
    plt.ylabel(column_y)
    plt.title(f'Scatter Plot cov {title}  - N={len(df)}') 

    plt.grid(True)  # Add grid (optional)
    plt.tight_layout()  # Adjust layout (optional)
    plt.show()

def simple_contour_plot(df, column_x, column_y, title):

    plt.figure(figsize=(8, 6))  # Define the size of the plot

    # Scatter plot with transparency
    plt.scatter(df[column_x], df[column_y], alpha=0.5, label='Data Points')

    # KDE plot (Kernel Density Estimation) for density visualization
    sns.kdeplot(df[column_x], df[column_y], cmap='Blues', shade=True, shade_lowest=False, label='Density')

    # Set labels and title
    plt.xlabel(column_x)
    plt.ylabel(column_y)

    # Count of subjects
    count_subjects = (len(df[[column_x, column_y]].dropna()))
    plt.title(f'Simple scatter Plot, cov {title}. Subjects: {count_subjects}')

    plt.legend()  # Show legend
    plt.grid(True)  # Add grid
    plt.tight_layout()
    plt.show()


def plot_with_polynomial_trendlines(df, column_x, column_y):

    plt.figure(figsize=(8, 6))  # Define the size of the plot

    # Scatter plot
    plt.scatter(df[column_x], df[column_y], alpha=0.5, label='Data Points')

    # Polynomial trendlines
    degrees = [1, 2, 3, 4, 5]  # Grados polinomiales a probar
    colors = ['green', 'blue', 'orange', 'purple', 'red']  # Colores para las líneas de tendencia

    for i, degree in enumerate(degrees):
        coefficients = np.polyfit(df[column_x], df[column_y], degree)
        polynomial_trendline = np.poly1d(coefficients)
        plt.plot(df[column_x], polynomial_trendline(df[column_x]), label=f'Degree {degree} Polynomial', color=colors[i])

    # Set labels and title
    plt.xlabel(column_x)
    plt.ylabel(column_y)
    plt.title('Scatter Plot with Polynomial Trendlines')

    plt.legend()  # Show legend
    plt.grid(True)  # Add grid
    plt.tight_layout()
    plt.show()


def calculate_metrics(df, column_x, column_y, degrees):

    for degree in degrees:
        # Ajustar modelo polinomial
        coefficients = np.polyfit(df[column_x], df[column_y], degree)
        polynomial_trendline = np.poly1d(coefficients)
        
        # Calcular predicciones del modelo
        y_pred = polynomial_trendline(df[column_x])
        
        # Calcular métricas
        mse = mean_squared_error(df[column_y], y_pred)
        r2 = r2_score(df[column_y], y_pred)
        
        # Mostrar métricas
        print(f"Degree {degree} Polynomial:")
        print(f"  MSE: {mse}")
        print(f"  R² : {r2}")
        print("-------------------------")

##############

def df_pair_regression(df, name_x, name_y, errors_plot=True):
    # Clean x, y from df
    x, y = x_y_drop_nan_df_pairs(df, name_x, name_y)
    # Compute main variables of the regression
    slope, intercept, r, pval, vstderr, fitted_y, errors = main_variables_regression(x,y)
    SSE, MSE, SST, R2, R2a = main_model_variables(y, errors)
    if errors_plot:
        errors_plots(y, fitted_y, errors)
    return slope, intercept, r, pval, vstderr, fitted_y, errors, SSE, MSE, SST, R2, R2a

def x_y_drop_nan_df_pairs(df, name_x, name_y):
    # remove nans from the pair:
    df_cleaned = df.dropna(subset=[name_x, name_y])
    x = df_cleaned[name_x] # covariable
    y = df_cleaned[name_y] # response variable
    return x, y

def main_variables_regression(x,y):
    slope , intercept , r , pval , stderr = stats.linregress(x,y)
    fitted_y = (x*slope+intercept)
    errors = y - fitted_y
    return slope, intercept , r, pval , stderr, fitted_y, errors


def main_model_variables(y, errors):
    SSE = sum( errors**2 )
    MSE = SSE/(len(y)-2)
    SST = sum( ( y - np.mean(y) )**2 )               # sum of square total
    R2 =  1 - SSE/SST                                # coefficient of determination
    R2a = 1 - ( (len(y)-1)/(len(y) - 2) ) * (1 - R2) # here, we have 2 parameters : intercept and slope
    print( """model 
    - SSE: {:.3f}
    - MSE: {:.3f}
    - R2:  {:.5f}
    - R2a: {:.5f}""".format(SSE,MSE,R2,R2a))
    return SSE, MSE, SST, R2, R2a

def errors_plots(y, fitted_y, errors):
    print('To ilustrate Strict exogeneity and Spertical error:')
    fig, ax = plt.subplots(ncols=3, figsize = (14,5))
    plot_errors_fitted_values(ax, fitted_y, errors)
    plot_count_errors(ax, errors)
    plot_real_fitted_values(ax, y, fitted_y)

def plot_errors_fitted_values(ax, fitted_y, errors):
    sns.scatterplot(x=fitted_y,y=errors, ax=ax[0], marker=".") # plot the errors along x
    ax[0].axhline(0,color='red')
    ax[0].set_ylabel('errors')
    ax[0].set_xlabel('fitted values')

def plot_count_errors(ax, errors):
    sns.histplot(errors, ax=ax[1]) # plot an histogram of the errors
    median_value = np.median(errors)
    max_counts = np.histogram(errors, bins='auto')[0].max()  # Maximum count of the histogram
    
    ax[1].vlines(median_value, 0, max_counts, color="red", linewidth=2)
    ax[1].set_xlabel('errors')
    # Add legend with additional information
    title_info = errors_basic_info(errors)
    ax[1].set_title(label=title_info)

def plot_real_fitted_values(ax, y, fitted_y):
    sns.scatterplot(x=fitted_y, y=y, ax=ax[2], marker=".") 
    ax[2].set_ylabel('real values')
    ax[2].set_xlabel('fitted values')

def errors_basic_info(errors):
    return [
        "error: mean {:.3f}".format(errors.mean()),
        "median {:.3f}".format(np.median(errors)), 
        "Q1 {:.3f}, Q3 {:.3f}".format(*np.quantile(errors, [0.25, 0.75]))
        ]

def scatter_plot_secondary_figs(df):

    list_param_a = [
        "tortuosity",
        "A tortuosity",
    ]  # , 'std diameter', 'V std diameter']
    list_param_b = [
        "Asc aorta min area"

    ]

    title = "no"
    multiple_scatterplots(df, list_param_a, list_param_b, title)

    ## Individuals:
    simple_scatter_plot(df, "A vascular density", "Asc aorta min area", title)

    simple_contour_plot(df, "A vascular density", "Asc aorta min area", title)

    df_cleaned = df.dropna(subset=["A vascular density", "Asc aorta min area"])
    plot_with_trendlines(
        df_cleaned, "A vascular density", "Asc aorta min area", title
    )


    # Plot polinomic tendency lines
    plot_with_polynomial_trendlines(
        df_cleaned, "A vascular density", "Asc aorta min area"
    )

    # Polinomic degress to try
    degrees_to_test = [1, 2, 3, 4, 5]

    calculate_metrics(
        df_cleaned, "A vascular density", "Asc aorta min area", degrees_to_test
    )

    ###  El R² es bastante bajo para todos los grados => los modelos polinomiales no están explicando bien la variación en los datos. (Un R² cercano a 1 indicaría que el modelo explica bien la variación de los datos)
    # El MSE, medida del error cuadrático medio entre las predicciones y los valores reales. Ninguno de los grados polinomiales produce un ajuste preciso, ya que el MSE es alto.

    # # Linear assumptions
    ### To change the model used for the LR to adapt better with the current code
    errors_analysis = True
    if errors_analysis:
        name_x = "A vascular density"
        name_y = "Asc aorta min area"

        (
            slope,
            intercept,
            r,
            pval,
            vstderr,
            fitted_y,
            errors,
            SSE,
            MSE,
            SST,
            R2,
            R2a,
        ) = df_pair_regression(df, name_x, name_y, errors_plot=True)
