#################################################################
#################################################################
###############  ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, glob
import pandas as pd
import rpy2.robjects as robjects

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
sys.path.append('../scripts')
import Support as S
import PipelineMcf10aP100 as P
from cmapPy.pandasGEXpress import parse

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
rawDataFiles = glob.glob('../rawdata/p100/LINCS_P100_*.gct')

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-mcf10a-p100.R'
r = robjects.r
r.source(rSource)

#######################################################
#######################################################
########## S1. Expression Data
#######################################################
#######################################################

#############################################
########## 1. Prepare Expression Data
#############################################

@follows(mkdir('s1-expression_data.dir'))

@merge(rawDataFiles,
   	   's1-expression_data.dir/p100-data.txt')

def getData(infiles, outfile):

	# Create empty list
	dataframes = []

	# Loop through infiles
	for infile in infiles:

		# Read GCT
		gctFile = parse(infile)

		# Get dataframe
		p100_dataframe = gctFile.data_df

		# Get names
		row_metadata_dataframe = gctFile.row_metadata_df
		row_metadata_dataframe['label'] = row_metadata_dataframe['pr_gene_symbol'] + '_' + row_metadata_dataframe['pr_p100_phosphosite']
		row_label_dict = row_metadata_dataframe.to_dict()['label']

		# Fix names
		p100_dataframe.index = [row_label_dict[x] for x in p100_dataframe.index]
		p100_dataframe.index.name = 'phospho_symbol'

		# Append
		dataframes.append(p100_dataframe)

	# Merge
	merged_dataframe = pd.concat(dataframes, axis=1)

	# Save
	merged_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S2. Sample Metadata
#######################################################
#######################################################

#############################################
########## 1. Get Sample Metadata
#############################################

@follows(mkdir('s2-annotations.dir'))

@merge(rawDataFiles,
	   's2-annotations.dir/p100-sample_annotations.txt')

def getSampleAnnotations(infiles, outfile):

	# Get dataframe
	sample_annotation_dataframe = pd.concat([parse(x).col_metadata_df for x in infiles])

	# Write
	sample_annotation_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')
