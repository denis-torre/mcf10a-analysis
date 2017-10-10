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
signatureMatrixfile = 's3-differential_expression.dir/l1000-signature_matrix.txt'

# Drugs
with open(drugFile) as openfile:
	drugs_dict = json.loads(openfile.read())
	drugs = drugs_dict['drugs']
	selected_drugs = drugs_dict['selected_drugs']

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
########## S6. Enrichment
#######################################################
#######################################################

#############################################
########## 1. Run Enrichr
#############################################

@follows(mkdir('s5-enrichr.dir'))

@files(mergeDifferentialExpression,
	   's5-enrichr.dir/l1000-enrichr_list_ids.txt')

def runEnrichr(infile, outfile):

	# Get signature dataframe
	signature_dataframe = pd.read_table(infile, index_col='gene_symbol')

	# Define dataframe list
	dataframes = []

	# Upload genesets
	for signature in signature_dataframe.columns:
		enrichr_id_dataframe = pd.DataFrame(S.uploadToEnrichr(signature_dataframe, signature, 100))
		enrichr_id_dataframe['signature'] = signature
		dataframes.append(enrichr_id_dataframe)

	# Concatenate
	merged_dataframe = pd.concat(dataframes)

	# Write
	merged_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Get results
#############################################

@transform(runEnrichr,
		   suffix('list_ids.txt'),
		   'results.txt')

def getEnrichmentResults(infile, outfile):

	# Get id dataframe
	enrichr_id_dataframe = pd.read_table(infile)

	# Define dataframe list
	dataframes = []

	# Loop through genesets
	for index, rowData in enrichr_id_dataframe.iterrows():

		# Get enrichment results
		enrichment_results_dataframe = S.getEnrichmentResults(rowData['userListId'], gene_set_library='GO_Biological_Process_2015')
		enrichment_results_dataframe['signature'] = rowData['signature']
		enrichment_results_dataframe['geneset'] = rowData['geneset']

		# Append
		dataframes.append(enrichment_results_dataframe)

	# Concatenate
	merged_dataframe = pd.concat(dataframes)

	# Write
	merged_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 3. Cast
#############################################

@follows(mkdir('s5-enrichr.dir/tables'))

@subdivide(getEnrichmentResults,
	   	   formatter(),
	   	   's5-enrichr.dir/tables/l1000-*-table.txt',
	   	   's5-enrichr.dir/tables/l1000-')

def castEnrichmentResults(infile, outfiles, outfileRoot):

	# Get id dataframe
	enrichment_dataframe = pd.read_table(infile)

	# Get logP
	enrichment_dataframe['logFDR'] = -np.log10(enrichment_dataframe['FDR'])

	# Get subset
	geneset_enrichment_dataframes = {geneset: pd.pivot_table(enrichment_dataframe.query('geneset == "{geneset}"'.format(**locals())), index='term_name', columns='signature', values='logFDR').fillna(0) for geneset in enrichment_dataframe['geneset'].unique()}

	# Save
	for key, value in geneset_enrichment_dataframes.iteritems():

		# Get outfile
		outfile = '{outfileRoot}{key}-table.txt'.format(**locals())

		# Save
		value.to_csv(outfile, sep='\t')

#############################################
########## 4. Merge
#############################################

@merge(castEnrichmentResults,
	   's5-enrichr.dir/l1000-enrichment_table.txt')

def mergeEnrichmentResults(infiles, outfile):

	# Read data
	enrichment_tables = {os.path.basename(x).split('-')[1]: pd.read_table(x, index_col='term_name') for x in infiles}

	# Filter terms
	filtered_terms = {key: value.apply(np.sum, 1).sort_values(ascending=False).index for key, value in enrichment_tables.iteritems()}

	# Filter data
	filtered_enrichment_tables = {key: value.loc[filtered_terms[key]] for key, value in enrichment_tables.iteritems()}

	# Change sign
	filtered_enrichment_tables['downregulated'] = -filtered_enrichment_tables['downregulated']

	# Get terms
	terms = set(filtered_enrichment_tables['upregulated'].index).intersection(set(filtered_enrichment_tables['downregulated'].index))

	# Get data dict
	enrichment_dict = {x:{} for x in filtered_enrichment_tables['upregulated'].index}

	# Fill
	for term in terms:
		for signature in filtered_enrichment_tables['upregulated'].columns:
			direction = 'upregulated' if filtered_enrichment_tables['upregulated'].loc[term, signature] >= abs(filtered_enrichment_tables['downregulated'].loc[term, signature]) else 'downregulated'
			enrichment_dict[term][signature] = filtered_enrichment_tables[direction].loc[term, signature]

	# Convert to dataframe
	enrichment_dataframe = pd.DataFrame(enrichment_dict).fillna(0).T
	enrichment_dataframe.index.name = 'term_name'

	# Write
	enrichment_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S6. Browser Data
#######################################################
#######################################################

#############################################
########## 1. Differential Expression
#############################################

def browserJobs():
	infile = signatureMatrixfile
	for drug in selected_drugs:
		outfile = 's6-browser_data.dir/{drug}/l1000-{drug}-browser_genes.json'.format(**locals())
		if not os.path.exists(os.path.dirname(outfile)):
			os.makedirs(os.path.dirname(outfile))
		yield [infile, outfile]

@follows(mkdir('s6-browser_data.dir'))

@files(browserJobs)

def getBrowserGenes(infile, outfile):

	# Read data
	signature_dataframe = pd.read_table(infile, index_col='gene_symbol')

	# Get drug
	drug = os.path.basename(outfile).split('-')[1]

	# Filter dataframe
	filtered_signature_dataframe = signature_dataframe[[x for x in signature_dataframe.columns if drug in x]]

	# Get top genes
	top_genes = abs(filtered_signature_dataframe).apply(np.sum, 1).sort_values(ascending=False).index.tolist()[:25]

	# Write
	with open(outfile, 'w') as openfile:
		openfile.write(json.dumps(top_genes))

#############################################
########## 2. Get Browser Signatures
#############################################

@transform(getBrowserGenes,
		   suffix('genes.json'),
		   add_inputs(mergeDifferentialExpression),
		   'signatures.txt')

def getBrowserSignatures(infiles, outfile):

	# Split infiles
	genesFile, signatureFile = infiles

	# Read genes
	with open(genesFile) as openfile:
		genes = json.loads(openfile.read())

	# Get drug
	drug = os.path.basename(outfile).split('-')[1]

	# Read data
	signature_dataframe = pd.read_table(signatureFile).set_index('gene_symbol', drop=False)

	# Filter
	filtered_signature_dataframe = signature_dataframe.loc[genes, ['gene_symbol']+[x for x in signature_dataframe.columns if drug in x]]

	# Melt
	melted_signature_dataframe = pd.melt(filtered_signature_dataframe, id_vars='gene_symbol', var_name='signature')

	# Split data
	signature_data = [x.split('_') for x in melted_signature_dataframe['signature']]
	melted_signature_dataframe['drug'] = [x[0] for x in signature_data]
	melted_signature_dataframe['concentration'] = [float(x[1]) for x in signature_data]
	melted_signature_dataframe['timepoint'] = [int(x[2]) for x in signature_data]

	# Subset
	filtered_melted_signature_dataframe = melted_signature_dataframe.query('drug=="{drug}"'.format(**locals())).sort_values(['gene_symbol', 'timepoint', 'concentration']).query('timepoint == 24')

	# Save
	filtered_melted_signature_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 3. Get Browser Coexpression
#############################################

@follows(getBrowserSignatures)

@transform(getBrowserGenes,
		   suffix('genes.json'),
		   add_inputs(mergeGeneCoexpression),
		   'coexpression.txt')

def getBrowserCoexpression(infiles, outfile):

	# Split infiles
	genesFile, coexpressionFile = infiles

	# Read genes
	with open(genesFile) as openfile:
		genes = json.loads(openfile.read())

	# Get drug
	drug = os.path.basename(outfile).split('-')[1]

	# Read data
	coexpression_dataframe = pd.read_table(coexpressionFile).set_index('source_gene_symbol', drop=False)

	# Filter
	genes_str = '("'+'", "'.join(genes)+'")'
	filtered_coexpression_dataframe = coexpression_dataframe[['source_gene_symbol', 'target_gene_symbol', drug]].rename(columns={drug: 'value'}).query('source_gene_symbol in {genes_str} and target_gene_symbol in {genes_str}'.format(**locals()))

	# Convert to json
	coexpression_data = [{'source': {'id': rowData['source_gene_symbol'], 'start': 0, 'end': 1}, 'target': {'id': rowData['target_gene_symbol'], 'start': 0, 'end': 1}, 'value': rowData['value']} for index, rowData in filtered_coexpression_dataframe.iterrows()]

	# Write
	filtered_coexpression_dataframe.to_csv(outfile, sep='\t')

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
