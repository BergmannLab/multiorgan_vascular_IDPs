##### MAIN RETINA PHENOTYPES

HOMOLOGOUS_LABELS_RED='DF_all,DF_artery,DF_vein,VD_orig_all,VD_orig_artery,VD_orig_vein,eq_CRAE,eq_CRVE,medianDiameter_all,medianDiameter_artery,medianDiameter_vein'
HOMOLOGOUS_NAMES_RED='tortuosity,A tortuosity,V tortuosity,vascular density,A vascular density,V vascular density,A central retinal eq,V central retinal eq,median diameter,A median diameter,V median diameter'

list_retina_homologous_red = list(HOMOLOGOUS_LABELS_RED.split(","))
list_retina_homologous_red_new = list(HOMOLOGOUS_NAMES_RED.split(","))

## OTHER OPTIONS
MAIN_LABELS='mean_angle_taa,mean_angle_tva,tau1_vein,tau1_artery,ratio_AV_DF,eq_CRAE,ratio_CRAE_CRVE,D_A_std,D_V_std,eq_CRVE,ratio_VD,VD_orig_artery,bifurcations,VD_orig_vein,medianDiameter_artery,medianDiameter_vein,ratio_AV_medianDiameter'
MAIN_NAMES='A temporal angle,V temporal angle,V tortuosity,A tortuosity,ratio tortuosity,A central retinal eq,ratio central retinal eq,A std diameter,V std diameter,V central retinal eq,ratio vascular density,A vascular density,bifurcations,V vascular density,A median diameter,V median diameter,ratio median diameter'

list_retina_IDPs = list(MAIN_LABELS.split(","))
list_retina_IDPs_new = list(MAIN_NAMES.split(","))

HOMOLOGOUS_LABELS_ALL='DF_all,DF_artery,DF_vein,VD_orig_all,VD_orig_artery,VD_orig_vein,CRAE,CRVE,eq_CRAE,eq_CRVE,median_CRAE,median_CRVE,medianDiameter_all,medianDiameter_artery,medianDiameter_vein,D_A_std,D_V_std,D_std,mean_intensity,std_intensity,nVessels'
HOMOLOGOUS_NAMES_ALL='tortuosity,A tortuosity,V tortuosity,vascular density,A vascular density,V vascular density,A central retinal eq2,V central retinal eq2,A central retinal eq,V central retinal eq,A central retinal median,V central retinal median,median diameter,A median diameter,V median diameter,A std diameter,V std diameter,std diameter,mean intensity,std intensity,num vessels'

list_retina_homologous = list(HOMOLOGOUS_LABELS_ALL.split(","))
list_retina_homologous_new = list(HOMOLOGOUS_NAMES_ALL.split(","))
