covariants_common= {
    '21003': 'Age', # Age Attended. it is more accurate
    #'21022': 'Age at recruitment',
    '31': 'Sex',
    '54': 'UK Biobank assessment centre',
    '50': 'Standing height'
    }

covariants_GPC = {
    '22009-0.1': 'PC1',
    '22009-0.2': 'PC2',
    '22009-0.3': 'PC3',
    '22009-0.4': 'PC4',
    '22009-0.5': 'PC5',
    '22009-0.6': 'PC6',
    '22009-0.7': 'PC7',
    '22009-0.8': 'PC8',
    '22009-0.9': 'PC9',
    '22009-0.10': 'PC10',
    '22009-0.11': 'PC11',
    '22009-0.12': 'PC12',
    '22009-0.13': 'PC13',
    '22009-0.14': 'PC14',
    '22009-0.15': 'PC15',
    '22009-0.16': 'PC16',
    '22009-0.17': 'PC17',
    '22009-0.18': 'PC18',
    '22009-0.19': 'PC19',
    '22009-0.20': 'PC20'
}

covariants_BP = { #not as main covar
    '4079': 'DBP', 
    '4080': 'SBP'
}

dict_all_covar = {}
dict_all_covar.update(covariants_common)
dict_all_covar.update(covariants_GPC)
dict_all_covar.update(covariants_BP)

list_all_covar_keys = list(dict_all_covar.keys())
list_all_covar_values = list(dict_all_covar.values())

# I will like to add diet (not sure if possible to do it accurately=> avoid), 
# sleep, exerice (no, but screens use, phone only small ss so not included), sun exposure (same=> avoid) if the data is relaible
# I do not include educational level since eye shape (reading effects) is already included


#####################################
#other_cov= {
    #'21002': 'Weight',
    #'21001': 'BMI',
    # '48': 'Waist circumference',  
    # '49': 'Hip circumference',
    # '22000': 'Genotype measurement batch', #no
    # '20161': 'Pack years of smoking', #no
    # '1558': 'Alcohol intake frequency', #no
    # '28693': 'Currently suffering from unrestful sleep', #yes, no, etc
    # '1200': 'Sleeplessness / insomnia', # categorial: rarely, always, etc
    # '1160': 'Sleep duration', #no
    # '1070':'Time spent watching television (TV)', #no
    # '1080':'Time spent using computer', #no
#}
#### Carotid measures:
# covariants_carotid={
#     '22682': 'Quality control indicator for IMT at 120 degrees',
#     '22683': 'Quality control indicator for IMT at 150 degrees',
#     '22684': 'Quality control indicator for IMT at 210 degrees',
#     '22685': 'Quality control indicator for IMT at 240 degrees'
# }

#### Brain measures: (I just add the ones from the anex from the paper: 'Confound modelling in UK Biobank brain imaging, Fidel Alfaro-Almagro, et al.')
## not such a paper in other organs. Confound familiy (Subject, adquisition, motion, table, non linear, crossed terms, time)
#covariants_brain={
    #'25000': 'Head size'#,
    ##'54': 'Acquisition site', (already included)
    ##'TBD': 'Batch', #TBD=to be decided
    ##'TBD': 'CMRR',
    ##'TBD': 'Protocol',
    ##'TBD': 'Service pack',
    ##'TBD': 'Scan ramp',
    ##'TBD': 'Scan cold head',
    ##'TBD': 'Scan head coil',
    ##'TBD': 'Scan misc',
    ##'TBD': 'Flipped SWI',
    #'26500': 'FS T2',
    #'25921': 'New eddy',
    #'25925': 'Scaling T1',
    #'25926': 'Scaling T2 FLAIR',
    #'25927': 'Scaling SWI',
    #'25928': 'Scaling dMRI',
    #'25929': 'Scaling rfMRI',
    #'25930': 'Scaling tfMRI',
    #'25923': 'TE rfMRI',
    #'25924': 'TE tfMRI',
    ##'TBD': 'Motion struct motion',
    #'18': 'DVARS',
    #'19': 'Head motion',
    #'10': 'Head motion ST',
    #'25756': 'Scan position X',
    #'25757': 'Scan position Y',
    #'25758': 'Scan position Z',
    #'25759': 'Scan table pos.',
    #'4': 'Eddy QC',
    ##'N/A': 'Non linear', # N/A= not applicable
    ##'N/A': 'Crossed terms',
    ##'53': 'Time acq time', #to split
    ##'53': 'Acq date' #to split
    #'53': 'Time and date acq' #to split
#}


#### Heart measures: Cardiac index and output seems more like output
# something else?
#covariants_heart={
    #'54': 'Acquisition site',
    #'22426':'Average heart rate', # number of heartbeats per minute at rest or during a specific period of time. beats per minute (bpm) 
    #'22427': 'Body surface area' # extent of a person's skin surface (m^2)
    #'22425': 'Cardiac index', # measure of cardiac function that takes into account cardiac output (the amount of blood the heart pumps per minute) and the patient's body surface area. It is used to assess the efficiency of the heart in relation to body size (L/min/m²).
    #'22424': 'Cardiac output' # Cardiac output is the total amount of blood the heart pumps per minute and is a fundamental indicator of cardiac function (L/min) 
#            }

#### Eye measures: 
# 1306	Refractometry	1+155
# 100099	Eye surgery	6
# 100016	Retinal optical coherence tomography	15+104
# 100017	Visual acuity	32
# 100015	Intraocular pressure	18
# covariants_eye = {
#     '5089': 'Astigmatism angle (left)',
#     '5088': 'Astigmatism angle (right)',
#     '5086': 'Cylindrical power (left)',
#     '5087': 'Cylindrical power (right)',
#     '5085': 'Spherical power (left)',
#     '5084': 'Spherical power (right)',
#     '5274': 'Vertex distance (left)',
#     '5215': 'Vertex distance (right)',
#     '7539': 'Pupil size (left)',
#     '7544': 'Pupil size (right)',
#     '5181': 'Ever had eye surgery',
#     '5187': 'Visual acuity measured (left)',
#     '5185': 'Visual acuity measured (right)',
#     '5265': 'Corneal resistance factor (left)',
#     '5257': 'Corneal resistance factor (right)',
#     '5262': 'Intra-ocular pressure, corneal-compensated (left)',
#     '5254': 'Intra-ocular pressure, corneal-compensated (right)',
#     '5263': 'Intra-ocular pressure, Goldmann-correlated (left)',
#     '5255': 'Intra-ocular pressure, Goldmann-correlated (right)'
# }


