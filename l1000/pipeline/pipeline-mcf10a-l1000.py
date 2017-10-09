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
import sys, json, os
import pandas as pd
import numpy as np
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
sys.path.append('../scripts')
import Support as S
import PipelineMcf10aL1000 as P
from cmapPy.pandasGEXpress import parse
from mcf10a import *

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
# Files
rawDataFile = '../rawdata/l1000/GSE70138_Broad_LINCS_Level2_GEX_n345976x978_2017-03-06.gctx'
geneAnnotationFile = '../rawdata/l1000/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt'
sampleAnnotationFile = '../rawdata/l1000/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt'
expressionDataFile = 's1-expression_data.dir/l1000-data.txt'
processedSampleAnnotationFile = 's2-annotations.dir/l1000-sample_annotations.txt'
drugFile = '../drugs.json'

# Drugs
with open(drugFile) as openfile:
	drugs = json.loads(openfile.read())['drugs']

##### 2. R Connection #####
rSource = 'pipeline/scripts/pipeline-mcf10a-l1000.R'
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
	sample_annotation_dataframe['pert_iname'] = [x.lower() for x in sample_annotation_dataframe['pert_iname']]
	sample_annotation_dataframe['pert_time'] = sample_annotation_dataframe['pert_time'].astype('int')
	sample_annotation_dataframe.index.name = 'sample_id'

	# Write
	sample_annotation_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S3. Differential Expression
#######################################################
#######################################################

#############################################
########## 1. Run Characteristic Direction
#############################################

def differentialExpressionJobs():
	group_count_dataframe = pd.read_table(processedSampleAnnotationFile).groupby(['pert_iname', 'pert_dose', 'pert_time']).size().rename('count').to_frame().reset_index()
	filtered_count_dataframe = group_count_dataframe.query('pert_iname != "dmso" and count > 2 and pert_iname in ("'+'", "'.join(drugs)+'")')
	infiles = [expressionDataFile, processedSampleAnnotationFile]
	for index, rowData in filtered_count_dataframe.iterrows():
		outfile = 's3-differential_expression.dir/cd/l1000-{pert_iname}_{pert_dose}_{pert_time}-differential_expression.txt'.format(**rowData)
		yield [infiles, outfile]

@follows(mkdir('s3-differential_expression.dir/cd'))

@files(differentialExpressionJobs)

def runDifferentialExpression(infiles, outfile):

	# Get infiles
	expressionDataFile, sampleAnnotationFile = infiles

	# Get expression data
	l1000_dataframe = pd.read_table(expressionDataFile, index_col='gene_symbol')

	# Get design
	drug, dose, timepoint = os.path.basename(outfile).split('-')[1].split('_')

	# Get annotation dataframe
	sample_annotation_dataframe = pd.read_table(sampleAnnotationFile, index_col='sample_id')
	sample_annotation_dataframe.index = [x.replace(':', '.') for x in sample_annotation_dataframe.index]

	# Get samples
	treated_samples = sample_annotation_dataframe.query('pert_iname == "{drug}" and pert_dose == {dose} and pert_time == {timepoint}'.format(**locals())).index.tolist()
	control_samples = sample_annotation_dataframe.query('pert_iname == "dmso" and pert_time == {timepoint}'.format(**locals())).index.tolist()

	# Run characteristic direction
	print 'Running {outfile}...'.format(**locals())
	cd_dataframe = pandas2ri.ri2py(r.runCharacteristicDirection(pandas2ri.py2ri(l1000_dataframe), treated_samples, control_samples))
	cd_dataframe.index.name = 'gene_symbol'

	# Write
	cd_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 2. Merge Results
#############################################

@merge(runDifferentialExpression,
	   's3-differential_expression.dir/l1000-signature_matrix.txt')

def mergeDifferentialExpression(infiles, outfile):

	# Define dataframe list
	dataframes = []

	# Loop through infiles
	for infile in infiles:

		# Get signature label
		signature = os.path.basename(infile).split('-')[1]

		# Read data
		dataframe = pd.read_table(infile, index_col='gene_symbol').rename(columns={'CD': signature})

		# Append
		dataframes.append(dataframe)

	# Concatenate
	# merged_dataframe = pd.concat(dataframes, axis=1)
	merged_dataframe = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True), dataframes)

	# Save
	merged_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S4. Coexpression
#######################################################
#######################################################

#############################################
########## 1. Gene coexpression by drug
#############################################

def coexpressionJobs():
	group_count_dataframe = pd.read_table(processedSampleAnnotationFile).groupby(['pert_iname']).size().rename('count').to_frame().reset_index()
	filtered_count_dataframe = group_count_dataframe.query('pert_iname != "dmso" and count > 2 and pert_iname in ("'+'", "'.join(drugs)+'")')
	infiles = [expressionDataFile, processedSampleAnnotationFile]
	for index, rowData in filtered_count_dataframe.iterrows():
		outfile = 's4-coexpression.dir/corr/l1000-{pert_iname}-coexpression.txt'.format(**rowData)
		yield [infiles, outfile]

@follows(mkdir('s4-coexpression.dir/corr'))

@files(coexpressionJobs)

def getGeneCoexpression(infiles, outfile):

	# Report
	print 'Doing {outfile}...'.format(**locals())

	# Get infiles
	expressionDataFile, processedSampleAnnotationFile = infiles

	# Get expression data
	l1000_dataframe = pd.read_table(expressionDataFile, index_col='gene_symbol')

	# Get drug
	drug = os.path.basename(outfile).split('-')[1]

	# Get annotation dataframe
	sample_annotation_dataframe = pd.read_table(processedSampleAnnotationFile, index_col='sample_id')

	# Get samples
	samples = sample_annotation_dataframe.query('pert_iname == "{drug}"'.format(**locals())).index

	# Get normalized subset
	drug_expression_dataframe = l1000_dataframe[samples]/(l1000_dataframe[samples].apply(np.sum))

	# Get correlation
	correlation_dataframe = drug_expression_dataframe.T.corr(method='spearman')
	np.fill_diagonal(correlation_dataframe.values, np.nan)
	correlation_dataframe.index.name = 'source_gene_symbol'
	correlation_dataframe.columns.name = 'target_gene_symbol'

	# Melt
	melted_correlation_dataframe = pd.melt(correlation_dataframe.reset_index(), id_vars='source_gene_symbol').dropna()

	# Write
	melted_correlation_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Merge Coexpression
#############################################

@merge(getGeneCoexpression,
	   's4-coexpression.dir/l1000-coexpression_matrix.txt')

def mergeGeneCoexpression(infiles, outfile):

	# Define dataframe list
	dataframes = []

	# Loop through infiles
	for infile in infiles:

		# Get signature label
		signature = os.path.basename(infile).split('-')[1]

		# Read data
		dataframe = pd.read_table(infile, index_col=['source_gene_symbol', 'target_gene_symbol']).rename(columns={'value': signature})

		# Append
		dataframes.append(dataframe)

	# Concatenate
	# merged_dataframe = pd.concat(dataframes, axis=1)
	merged_dataframe = reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True), dataframes)

	# Save
	merged_dataframe.to_csv(outfile, sep='\t')


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
pipeline_run([sys.argv[-1]], multiprocess=4, verbose=1)
print('Done!')
