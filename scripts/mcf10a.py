#################################################################
#################################################################
###############  MCF10A Analysis Support Scripts ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
import sys, json
import pandas as pd
import numpy as np

#######################################################
#######################################################
########## S1. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Get Design
#############################################

def getDesignDataframe(sample_metadata_dataframe, drug_dict, timepoint=24):

	# Define sample dict
	dataframe_list = []

	# Loop through drugs
	for drug, drug_concentration_dict in drug_dict.iteritems():

		# Loop through drug concentrations
		for concentration_level, concentration_bounds in drug_concentration_dict.iteritems():

			# Add samples
			sample_dataframe = sample_metadata_dataframe.query('pert_time == {timepoint} and pert_iname == "{drug}"'.format(**locals())).query('pert_dose >= {min} and pert_dose <= {max}'.format(**concentration_bounds)).copy()
			sample_dataframe['concentration_level'] = concentration_level
			dataframe_list.append(sample_dataframe)

	# Concatenate
	result_dataframe = pd.concat(dataframe_list)

	# Fix column names
	rename_dict = {'pert_dose': 'dose', 'pert_iname': 'drug', 'pert_time': 'timepoint', 'concentration_level': 'level'}
	result_dataframe = result_dataframe.rename(columns=rename_dict)[rename_dict.values()]

	# Return
	return result_dataframe