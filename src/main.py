import sys
import typer
import logging

from settings import DIR_UTILS, PYTHON_LOGGING_LEVEL

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)

app = typer.Typer()

sys.path.append(DIR_UTILS)

from repository.idps_retina_cov.create_dataset import create_dataset_csv
from repository.preprocessing_pheno.preprocessing_idps_cov import preproces_idps_cov
from repository.preprocessing_pheno.idps_resid import calculate_resid_from_idps_csv
from repository.plotting_pheno.phenotypic_plots import plot_phenotypic_analysis
from repository.plotting_genet.preprocessing_gcorr import preprocesing_ldsr_output
from repository.plotting_genet.plot_gcorr import plot_ldsr_genetic_corr
from repository.plotting_genet.plot_shared_genes_pathways import plot_genes_pathways_intersection

logger.info("I) Creating the dataset\n")
create_dataset_csv()

logger.info("II) Phenotypic preprocessing (rename columns and filter IDPs by 2nd instance)\n")
preproces_idps_cov()

logger.info("III) Compute the idps residuals\n")
calculate_resid_from_idps_csv()

logger.info("IV) Phenotypic plotting. Also possible to use jupyternotebook\n")
plot_phenotypic_analysis()

logger.info("V) Organizing the genetic correlation (individual phenotype data must be computed already)\n")
preprocesing_ldsr_output()
logger.info("VI) Plotting the genetic correlation (LDSR). Also possible to use jupyternotebook\n")
plot_ldsr_genetic_corr()

logger.info("VI) Plotting the intersection of genes and pathways accros phenotypes pairs. Also possible to use jupyternotebook\n")
plot_genes_pathways_intersection()

if __name__ == "__main":
    app()