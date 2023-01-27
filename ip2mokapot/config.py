# File Inputs
SQTS_DESCRIPTION = 'Path to SQTs files'
FASTAS_DESCRIPTION = 'Path to FASTA file used to generate the provided SQT files'
SEARCH_XML_DESCRIPTION = 'The path to the search.xml file used to generate the input SQT files. Note: Including a search.xml file will override the following arguments: enzyme_regex, missed_cleavage, min_length, semi'
DTASELECT_PARAMS_DESCRIPTION='Path to dtaSelect.params file, Note: Including a DTASelect.params file will override the following arguments: protein_fdr, peptide_fdr, psm_fdr, and min_peptides'
OUT_DESCRIPTION = 'The output DTASelect-filter.txt path'

# DTASelect.params arguments
PROTEIN_FDR_DESCRIPTION = 'The Protein level FDR'
PEPTIDE_FDR_DESCRIPTION = 'The Peptide level FDR'
PSM_FDR_DESCRIPTION = 'The PSM level FDR'
MIN_PEPTIDES_DESCRIPTION = 'The minimum number of peptides required to identify a protein'

# Search.xml arguments
ENZYME_REGEX_DESCRIPTION = 'The regex string to determine enzyme sites (should match original search settings!)'
ENZYME_TERM_DESCRIPTION = 'The regex term (should match original search settings!)'
MISSED_CLEAVAGE_DESCRIPTION = 'The number of internal peptide missed cleavages to allow (should match original search settings!)'
MIN_LENGTH_DESCRIPTION = 'The minimum peptide length (should match original search settings)'
MAX_LENGTH_DESCRIPTION = 'The maximum peptide length (should match original search settings)'
SEMI_DESCRIPTION = 'A flag for specifying enzymatic specificity. True if both ends follow adhere to the enzyme_regex, False otherwise'

# FASTA arguments
DECOY_PREFIX_DESCRIPTION = 'The decoy prefix found int FASTA and SQT files. For IP2 use Reverse_'

# Mokapot arguments
XGBOOST_DESCRIPTION = 'Use Xbgoost as the model for mokapots semi-supervised training loop. If false mokapot will use its Percolator like model.'
TEST_FDR_DESCRIPTION = 'The test FDR to use during semi-supervized training loop'
FOLDS_DESCRIPTION = 'number of K-Folds to preform'
WORKERS_DESCRIPTION = 'The number of workers (threads) to use for semi-supervized training loop'
MAX_ITER_DESCRIPTION = 'The number of iterations to preform during the semi-supervized training loop'
TIMSCORE_DESCRIPTION = 'A flag to denote PaSER SQT files. Setting this to True will 1) use timscores as a feature in mokapot, and 2) output a PaSER compatible DTASelect-filter.txt file'
MAX_SEED_DESCRIPTION = 'The random seed to use, for reproducibility.'

# IP2Mokapot arguments
MASS_ALIGNMENT_DESCRIPTION = 'Align masses within each sqt file, and recalculate ppm. Alignment uses a 1D polyfit fitted to the mass ppm drift vs retention time for peptides >= 95th percentile xcorr'
MAX_MLINE_DESCRIPTION = 'The maximum number of m lines to use for input to mokapot.'

# Help msg
HELP_DESCRIPTION = """
    
    Welcome to MokaFilter! This tutorial will guide you through your first time using our app.
    
    To begin, upload one or more SQT files. These files contain the raw peptide search matches (PSMs) from the Prolucid search engine.
    Next, upload one or more FASTA files. These files(s) should be the ones used generate the SQT files
    Next, upload the search.xml file, which specifies the Prolucid search parameters.
    Next upload the dtaselect.params file, which contains the filtering parameters.
    Once all of the necessary files have been uploaded, click run and Ip2MokaPot will generate a filtered results list of proteins, peptides, and spectra at the provided false discovery rate (FDR).
    
    
    
    Fondrie W. E. & Noble W. S. mokapot: Fast and Flexible Semisupervised Learning for Peptide Detection. J Proteome Res (2021) doi: 10.1021/acs.jproteome.0c01010. PMID: 33596079.
    """