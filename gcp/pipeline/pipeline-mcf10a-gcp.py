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
import sys
import pandas as pd
import rpy2.robjects as robjects

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
sys.path.append('../scripts')
import Support as S
import PipelineMcf10aGcp as P
from cmapPy.pandasGEXpress import parse

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
rawDataFile = '../rawdata/gcp/LINCS_GCP_Plate35_annotated_minimized_2016-06-09_13-14-47.gct'

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-mcf10a-gcp.R'
r = robjects.r
r.source(rSource)

#######################################################
#######################################################
########## S1. Process Data
#######################################################
#######################################################

#############################################
########## 1. Process Dataset
#############################################

@follows(mkdir('s1-data.dir'))

@files(rawDataFile,
	   's1-data.dir/gcp-data.txt')

def processData(infile, outfile):

	# Parse
	gctFile = parse(infile)

	# Get data
	gcp_dataframe = gctFile.data_df

	# Rename 
	gcp_dataframe = gcp_dataframe.rename(index=gctFile.row_metadata_df.to_dict()['pr_gcp_histone_mark']).drop('H3NORM(41-49)')

	# Save
	gcp_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S2. Annotations
#######################################################
#######################################################

#############################################
########## 1. Sample Annotations
#############################################

@follows(mkdir('s2-annotations.dir'))

@files(rawDataFile,
	   's2-annotations.dir/gcp-sample_annotations.txt')

def getSampleAnnotations(infile, outfile):

	# Parse
	sample_metadata_dataframe = parse(infile).col_metadata_df

	# Write
	sample_metadata_dataframe.to_csv(outfile, sep='\t')

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
