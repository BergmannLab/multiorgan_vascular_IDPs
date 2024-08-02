# UKBB Vascular Image Derived Phenotypes (IDPs) Across Different Organs

## Project Description
This project aims to compare vascular Image Derived Phenotypes (IDPs) from different organs using UKBiobank data. It includes both phenotypic and genetic analyses:

- **Phenotypic Analysis**: Filters data by the first visit, computes the residuals of each IDP with respect to covariates, and plots the correlations of these residuals.
- **Genetic Analysis**: Uses individual phenotypic-pair correlations (previously computed using LDSR) to plot genetic correlations across multiple phenotypes. Additionally, it examines the number of significant genes and pathways shared among the phenotypes using individual phenotypic genes and pathways (previously computed using PascalX), with corresponding plots.

## Requirements
- Python version 3.8.0
- Required packages are listed in `requirements.txt`

## Configuration
Create a `.env` file in the `/VasculatureMultiOrgan` directory with the variables specified in the `example_conf_file.py`.

## Project Structure
``` 
├── src
│ ├── init.py
│ ├── main.sh
│ ├── main.py
│ ├── settings.py
│ ├── requirements.txt
│ ├── repository
│ │ ├── idps_retina_cov
│ │ ├── preprocessing_pheno
│ │ ├── plotting_pheno
│ │ ├── plotting_genet
│ ├── utils
│ │ ├── create_dataset
│ │ ├── data_information
│ │ ├── plotting
│ │ ├── preprocessing_pheno
│ │ ├── preprocessing_genet
│ ├── jupyternotebooks
├── requirements.txt
├── .env
├── README.md
├── .gitignore
``` 

## Documentation
1. Clone the repository.
2. Set up your `.env` file as described.
3. Ensure you are using Python 3.8.0 and have installed the required packages.
4. Create a `output/` folder in the `/VasculatureMultiOrgan/src/` directory.
5. Add your Python version to `main.sh` and run the script using `bash main.sh`.

### Output
The output will include:
- Phenotypic and genetic plots saved in the `output/` directory.
- CSV files with the phenotypic data saved in the directory specified in your `.env` file.

---