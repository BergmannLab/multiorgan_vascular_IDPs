####### BRAIN ######
# CBF is defined as the blood volume that flows per unit mass per unit time brain tissue and is typically expressed units of ml blood / ( 100 g tissue min )
dict_CBF = {
    "24380": "Mean CBF Caudate (L)",
    "24379": "Mean CBF Caudate (R)",
    "24363": "Mean CBF Cerebellum grey matter (L)",
    "24362": "Mean CBF Cerebellum grey matter (R)",
    "24365": "Mean CBF Frontal Lobe grey matter (L)",
    "24364": "Mean CBF Frontal Lobe grey matter (R)",
    "24373": "Mean CBF Internal Carotid Artery vascular territory grey matter (L)",
    "24372": "Mean CBF Internal Carotid Artery vascular territory grey matter (R)",
    "24367": "Mean CBF Occipital Lobe grey matter (L)",
    "24366": "Mean CBF Occipital Lobe grey matter (R)",
    "24369": "Mean CBF Parietal Lobe grey matter (L)",
    "24368": "Mean CBF Parietal Lobe grey matter (R)",
    "24382": "Mean CBF Putamen (L)",
    "24381": "Mean CBF Putamen (R)",
    "24371": "Mean CBF Temporal Lobe grey matter (L)",
    "24370": "Mean CBF Temporal Lobe grey matter (R)",
    "24384": "Mean CBF Thalamus (L)",
    "24383": "Mean CBF Thalamus (R)",
    "24374": "Mean CBF VertebroBasilar Arteries vascular territories grey matter",
    "24377": "Mean CBF cerebrum white matter (L)",
    "24378": "Mean CBF cerebrum white matter (R)",
    "24376": "Mean CBF cerebrum white matter and >50% cerebral partial volume",
    "24361": "Mean CBF cortex grey matter",
    "24360": "Mean CBF whole brain grey matter",
    "24375": "Mean CBF whole brain white matter"
    }

# Arterial transit time (ATT)
dict_ATT={
            "24405":	"Mean ATT Caudate (L)",
            "24404":	"Mean ATT Caudate (R)",
            "24388":	"Mean ATT Cerebellum grey matter (L)",
            "24387":	"Mean ATT Cerebellum grey matter (R)",
            "24390":	"Mean ATT Frontal Lobe grey matter (L)",
            "24389":	"Mean ATT Frontal Lobe grey matter (R)",
            "24398":	"Mean ATT Internal Carotid Artery vascular territory grey matter (L)",
            "24397":	"Mean ATT Internal Carotid Artery vascular territory grey matter (R)",
            "24392":	"Mean ATT Occipital Lobe grey matter (L)",
            "24391":	"Mean ATT Occipital Lobe grey matter (R)",
            "24394":	"Mean ATT Parietal Lobe grey matter (L)",
            "24393":	"Mean ATT Parietal Lobe grey matter (R)",
            "24407":	"Mean ATT Putamen (L)",
            "24406":	"Mean ATT Putamen (R)",
            "24396":	"Mean ATT Temporal Lobe grey matter (L)",
            "24395":	"Mean ATT Temporal Lobe grey matter (R)",
            "24409":	"Mean ATT Thalamus (L)",
            "24408":	"Mean ATT Thalamus (R)",
            "24399":	"Mean ATT VertebroBasilar Arteries vascular territories grey matter",
            "24402":	"Mean ATT cerebrum white matter (L)",
            "24403":	"Mean ATT cerebrum white matter (R)",
            "24401":	"Mean ATT cerebrum white matter and >50% cerebral partial volume",
            "24386":	"Mean ATT cortex grey matter",
            "24385":	"Mean ATT whole brain grey matter",
            "24400":	"Mean ATT whole brain white matter"}

dict_deletions_brain = {
            "24486":	"Deep WMH total volume",
            "24485":	"Peri-ventricular WMH",
            "25781":	"WMH total volume"}

dict_brain_vessel={
            "26550":"Brain mean vessel intensity (L)",
            "26581":"Brain mean vessel intensity (R)",
            "26566":"Brain vessel volume (L)",
            "26597":"Brain vessel volume (R)"}

dict_vascular_brain = dict_CBF.copy()
dict_vascular_brain.update(dict_ATT)
dict_vascular_brain.update(dict_deletions_brain)
dict_vascular_brain.update(dict_brain_vessel)

list_keys_brain = list(dict_vascular_brain.keys())
list_values_brain = list(dict_vascular_brain.values())


####### HEART ######
dict_heart_vessel={
            "24120":"Asc aorta distensibility",
            "24119":"Asc aorta min area",
            "24123":"Desc aorta distensibility",
            "24122":"Desc aorta min area",
            "24118":"Asc aorta max area",
            "24121":"Desc aorta max area"
            }

dict_heart_vessel_2={
            "31079":"Diastole mean area of Asc aorta",
            "31080":"Proximal Desc aorta Diastole mean area",
            "31083":"Diastole std of Asc aorta mean",
            "31084":"Proximal Desc aorta mean Diastole std",
            "31077":"Systole mean area of Asc aorta",
            "31078":"Proximal Desc aorta Systole mean area",
            "31081":"Asc aorta area Systole std",
            "31082":"Proximal Desc aorta area Systole std",
            "31128":"Asc aorta area mean abs deviation",
            "31129":"Proximal Desc aorta area mean abs deviation",
            }

dict_heart_functional={
            "24113":"LA ejection fraction",
            "22420":"LV ejection fraction2",
            "24103":"LV ejection fraction",
            "24117":"RA ejection fraction",
            "24109":"RV ejection fraction",
            "24112":"LA stroke volume",
            "24102":"LV stroke volume",
            "24116":"RA stroke volume",
            "24108":"RV stroke volume"
            }

dict_vascular_heart = dict_heart_vessel.copy()
dict_vascular_heart.update(dict_heart_vessel_2)
dict_vascular_heart.update(dict_heart_functional)

list_keys_heart = list(dict_vascular_heart.keys())
list_values_heart = list(dict_vascular_heart.values())

####### CAROTID ######
dict_IMT={
            "22672":"Max carotid IMT at 120 degrees",
            "22675":"Max carotid IMT at 150 degrees",
            "22678":"Max carotid IMT at 210 degrees",
            "22681":"Max carotid IMT at 240 degrees",
            "22671":"Mean carotid IMT at 120 degrees",
            "22674":"Mean carotid IMT at 150 degrees",
            "22677":"Mean carotid IMT at 210 degrees",
            "22680":"Mean carotid IMT at 240 degrees",
            "22670":"Min carotid IMT at 120 degrees",
            "22673":"Min carotid IMT at 150 degrees",
            "22676":"Min carotid IMT at 210 degrees",
            "22679":"Min carotid IMT at 240 degrees"
        }

dict_vascular_ultrasound = dict_IMT.copy()
#dict_vascular.update(dict_brain_vessel)

list_keys_carotid = list(dict_vascular_ultrasound.keys())
list_values_carotid = list(dict_vascular_ultrasound.values())
