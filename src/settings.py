import logging
import os

from dotenv import load_dotenv

load_dotenv("/SSD/home/sofia/VasculatureMultiOrgan/.env")
PYTHON_LOGGING_LEVEL = os.getenv("PYTHON_LOGGING_LEVEL", logging.INFO)

DIR_UTILS = os.getenv("DIR_UTILS")
DIR_UKBB = os.getenv("DIR_UKBB")
DIR_FILE_EID_MAP = os.getenv("DIR_FILE_EID_MAP")
DIR_FILE_RETINA_IDPS = os.getenv("DIR_FILE_RETINA_IDPS")
DIR_LDSR = os.getenv("DIR_LDSR")
DIR_OUTPUT = os.getenv('DIR_OUTPUT')
DIR_OTHER_GENES_PATH = os.getenv("DIR_OTHER_GENES_PATH")
DIR_GENES_RET = os.getenv("DIR_GENES_RET")
OLD_EID = os.getenv("OLD_EID")
NEW_EID = os.getenv("NEW_EID")
INTERM_FOLDER = os.getenv("INTERM_FOLDER")
DATE_USED = os.getenv("DATE_USED")
FILE_IDPS = os.getenv("FILE_IDPS")
FILE_GENETIC_INFO = os.getenv("FILE_GENETIC_INFO")
HIST_FILE = os.getenv("HIST_FILE")
END_RAW_COV = os.getenv("END_RAW_COV")
END_RAW_IDPS_COV_RETINA = os.getenv("END_RAW_IDPS_COV_RETINA")
IDP_RETINA_USED = os.getenv("IDP_RETINA_USED")
IDP_VASCULAR_USED = os.getenv("IDP_VASCULAR_USED")
COVARIATES = os.getenv('COVARIATES')
BP_USED = os.getenv('BP_USED')
TYPE_BP = os.getenv('TYPE_BP')
TYPE_OF_FIGS = os.getenv('TYPE_OF_FIGS')
SAVE_FIGURES = os.getenv('SAVE_FIGURES')
PLOT_SCATTER_FIGS = os.getenv('PLOT_SCATTER_FIGS')

Z_SCORE_APPLIED = os.getenv("Z_SCORE_APPLIED")
FILTER_OUTLIERS = os.getenv("FILTER_OUTLIERS")
NUM_STD = os.getenv("NUM_STD")
L_GEN_PATH = os.getenv("L_GEN_PATH")

# Figs names:
SQUARE_FIG = os.getenv("SQUARE_FIG")
RECT_FIG = os.getenv("RECT_FIG")
DOBLE_FIG = os.getenv("DOBLE_FIG")

# P values parameters:
P_VAL_GENES = os.getenv("P_VAL_GENES")
P_VAL_PATHS = os.getenv("P_VAL_PATHS")
# Multiple testing:
SOFT_MULTIPLE_TESTING = os.getenv("SOFT_MULTIPLE_TESTING")
MULTIPLE_TESTING = os.getenv("MULTIPLE_TESTING")
ALPHA_1 = float(os.getenv("ALPHA_1"))
ALPHA_2 = float(os.getenv("ALPHA_2"))

# Other parameters:
N_NAMES_SHOWN = os.getenv("N_NAMES_SHOWN")
PLOT_HISTOGRAMS = os.getenv("PLOT_HISTOGRAMS")
# Num genes/path names to show