import sys
import logging
import utils.data_information.data_IDPs_info as d_
from settings import PYTHON_LOGGING_LEVEL, DIR_UTILS

logger = logging.getLogger(__name__)
logging.basicConfig(level=PYTHON_LOGGING_LEVEL)

sys.path.append(DIR_UTILS)


dict_homologous = d_.dict_CBF.copy()
dict_homologous.update(d_.dict_ATT)
dict_homologous.update(d_.dict_deletions_brain)
dict_homologous.update(d_.dict_brain_vessel)
dict_homologous.update(d_.dict_IMT)
dict_homologous.update(d_.dict_heart_vessel)
dict_homologous.update(d_.dict_heart_vessel_2)
dict_homologous.update(d_.dict_vascular_heart)

dict_avg_brain = {
    'Mean CBF': 'Mean CBF', 
    'Mean ATT': 'Mean ATT'}

dict_avg_IMT = {
    'Min carotid IMT': 'Min carotid IMT',
    'Meann carotid IMT': 'Mean carotid IMT',
    'Max carotid IMT': 'Max carotid IMT'
}

dict_homologous_averag = dict_avg_brain.copy()
dict_homologous_averag.update(d_.dict_deletions_brain)
dict_homologous_averag.update(d_.dict_brain_vessel)
dict_homologous_averag.update(dict_avg_IMT)
dict_homologous_averag.update(d_.dict_vascular_heart)


## Manually filtered
dict_final_IDPs = {
            "25781":"WMH total volume",
            #"26550":"Brain mean vessel intensity (L)",
            #"26581":"Brain mean vessel intensity (R)",
            #"26566":"Brain vessel volume (L)",
            #"26597":"Brain vessel volume (R)",
            'Min carotid IMT':'Min carotid IMT',
            "24120":"Asc aorta distensibility",
            "24119":"Asc aorta min area",
            "24123":"Desc aorta distensibility",
            "24122":"Desc aorta min area" #,
           # "31083":"Diastole standard deviation of ascending aorta mean",
           # "31084":"Diastole standard deviation of proximal descending aorta mean",
           # "31081":"Systole standard deviation of ascending aorta area",
           # "31082":"Systole standard deviation of proximal descending aorta area",
           # "31128":"Mean absolute deviation of ascending aorta area",
           # "31129":"Mean absolute deviation of proximal descending aorta area",
        }