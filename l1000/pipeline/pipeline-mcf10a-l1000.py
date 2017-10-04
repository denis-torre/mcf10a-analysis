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
import numpy as np
import rpy2.robjects as robjects

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
sys.path.append('../scripts')
import Support as S
import PipelineMcf10aL1000 as P
from cmapPy.pandasGEXpress import parse

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
# Files
rawDataFile = '../rawdata/l1000/GSE70138_Broad_LINCS_Level2_GEX_n345976x978_2017-03-06.gctx'
geneAnnotationFile = '../rawdata/l1000/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt'
sampleAnnotationFile = '../rawdata/l1000/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt'

# Drugs
drugs = ['DMSO', 'trametinib', 'alpelisib', 'vorinostat', 'neratinib', 'paclitaxel', 'palbociclib', 'etoposide', 'dasatinib']

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-mcf10a-l1000.R'
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

@follows(mkdir('s1-expression_data.dir'))

@files([rawDataFile, geneAnnotationFile, sampleAnnotationFile],
	   's1-expression_data.dir/l1000-data.txt')

def getExpressionData(infiles, outfile):

	# Get infiles
	rawDataFile, geneAnnotationFil, sampleAnnotationFile = infiles

	# Get object
	gctxData = parse(rawDataFile)

	# Get expression dataframe
	expression_dataframe = gctxData.data_df

	# Get gene annotations
	gene_annotation_dict = pd.read_table(geneAnnotationFile).set_index('pr_gene_id').to_dict()['pr_gene_symbol']

	# Fix gene symbols
	expression_dataframe.index = [gene_annotation_dict[x] for x in expression_dataframe.index]

	# Get sample annotations
	mcf10a_samples = pd.read_table(sampleAnnotationFile).query('cell_id=="MCF10A"').query('pert_iname in ("'+'", "'.join(drugs)+'")')['inst_id']

	# Get subset
	expression_dataframe_subset = expression_dataframe[mcf10a_samples]

	# Add index
	expression_dataframe_subset.index.name = 'gene_symbol'

	# Save data
	expression_dataframe_subset.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S2. Annotations
#######################################################
#######################################################

#############################################
########## 1. Sample Annotations
#############################################

@follows(mkdir('s2-annotations.dir'))

@files([getExpressionData, sampleAnnotationFile],
	   's2-annotations.dir/l1000-sample_annotations.txt')

def getSampleAnnotations(infiles, outfile):

	# Get infiles
	expressionInfile, sampleAnnotationFile = infiles

	# Get MCF10A samples
	mcf10a_samples = pd.read_table(expressionInfile, index_col='gene_symbol').columns

	# Get annotation dataframe
	sample_annotation_dataframe = pd.read_table(sampleAnnotationFile, index_col='inst_id').loc[mcf10a_samples].replace(-666, np.nan).replace('-666', np.nan)
	sample_annotation_dataframe['pert_time'] = sample_annotation_dataframe['pert_time'].astype('int')

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
