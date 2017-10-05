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
import sys, json
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
from mcf10a import *

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
rawDataFile = '../rawdata/gcp/LINCS_GCP_Plate35_annotated_minimized_2016-06-09_13-14-47.gct'
processedDataFile = 's1-data.dir/gcp-data.txt'
sampleAnnotationFile = 's2-annotations.dir/gcp-sample_annotations.txt'
drugFile = '../drugs.json'

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
	sample_metadata_dataframe['pert_iname'] = [x.lower() for x in sample_metadata_dataframe['pert_iname']]

	# Write
	sample_metadata_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S3. Differential Expression
#######################################################
#######################################################

#############################################
########## 2. Get Design
#############################################

@follows(mkdir('s3-differential_expression.dir'))

@merge([getSampleAnnotations, drugFile],
	   's3-differential_expression.dir/gcp-design.txt')

def getDesign(infiles, outfile):

	# Get infiles
	sampleAnnotationFile, drugFile = infiles

	# Read sample annotations
	sample_metadata_dataframe = pd.read_table(sampleAnnotationFile, index_col='cid').fillna(0)#.rename(columns={''})

	# Read drug dict
	with open(drugFile) as openfile:
		drug_dict = json.loads(openfile.read())

	# Get design
	design_dataframe = getDesignDataframe(sample_metadata_dataframe, drug_dict)
	
	# Save
	design_dataframe.to_csv(outfile, sep='\t')

# #############################################
# ########## 1. Run Differential Expression
# #############################################

# def differentialExpressionJobs():
# 	drugs = pd.read_table(sampleAnnotationFile)['pert_iname'].unique()
# 	for drug in drugs:
# 		infiles = [processedDataFile, sampleAnnotationFile]
# 		outfile = 's3-differential_expression.dir/{drug}-gcp_differential_expression.txt'.format(**locals())
# 		yield [infiles, outfile]

# @follows(mkdir('s3-differential_expression.dir'))

# @files(differentialExpressionJobs)

# def runDifferentialExpression(infiles, outfile):

# # Get infiles
# processedDataFile, sampleAnnotationFile = infiles

# # Read data
# gcp_dataframe = pd.read_table(processedDataFile, index_col='rid')

# # Read sample annotations
# sample_metadata_dataframe = pd.read_table(sampleAnnotationFile)

# # Get drug
# drug = os.path.basename(outfile)[:-len('-gcp_differential_expression.txt')]

# # Get samples
# treated_samples = sample_metadata_dataframe.query('pert_time == 24 and pert_iname == "{drug}" and '.format(**locals()))['cid']

# print infiles, outfile



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
