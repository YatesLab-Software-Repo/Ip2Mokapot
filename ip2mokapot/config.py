# File Inputs
SQTS_DESCRIPTION = 'Path to .sqt files and/or folders containing .sqt files. These files represent the output from prolucid search engine and are used as input for further analysis.'
FASTAS_DESCRIPTION = 'Path to the .fasta file(s) and/or folder(s) containing .fasta files. These files contain the protein sequences used in the search engine.'
SEARCH_XML_DESCRIPTION = 'The path to the search.xml file used to generate the input SQT files. It contains information about the parameters used in the search engine. Note: including this file will override the following arguments: enzyme_regex, missed_cleavage, min_length, semi.'
DTASELECT_PARAMS_DESCRIPTION = 'Path to the dtaSelect.params file. This file contains the parameters used by dtaSelect to filter the results from the search engine. Note: including this file will override the following arguments: protein_fdr, peptide_fdr, psm_fdr, and min_peptides.'
OUT_DESCRIPTION = 'The output DTASelect-filter.txt path. This file will contain the filtered results from the search engine.'

# DTASelect.params arguments
PROTEIN_FDR_DESCRIPTION = 'The false discovery rate (FDR) for proteins.'
PEPTIDE_FDR_DESCRIPTION = 'The false discovery rate (FDR) for peptides.'
PSM_FDR_DESCRIPTION = 'The false discovery rate (FDR) for PSMs (peptide-spectrum matches).'
MIN_PEPTIDES_DESCRIPTION = 'The minimum sequence count required to identify a protein.'

# Search.xml arguments
ENZYME_REGEX_DESCRIPTION = 'The regular expression (regex) string used to determine the enzyme cleavage sites. This should match the settings used in the original search engine.'
ENZYME_TERM_DESCRIPTION = 'A flag indicating whether the enzyme cleaves at the C-terminus (True) or N-terminus (False) of the protein. This should match the settings used in the original search engine.'
MISSED_CLEAVAGE_DESCRIPTION = 'The number of internal peptide missed cleavages allowed in the search. This should match the settings used in the original search engine.'
MIN_LENGTH_DESCRIPTION = 'The minimum length of a peptide in the search. This should match the settings used in the original search engine.'
MAX_LENGTH_DESCRIPTION = 'The maximum length of a peptide in the search. This should match the settings used in the original search engine.'
SEMI_DESCRIPTION = 'A flag for specifying the enzymatic specificity of the search. True if both ends of a peptide follow the enzyme_regex, False otherwise. This should match the settings used in the original search engine.'

# FASTA arguments
DECOY_PREFIX_DESCRIPTION = 'The prefix used to identify decoy sequences in the FASTA and SQT files. For IP2, the decoy prefix is "Reverse_".'
# Mokapot arguments
XGBOOST_DESCRIPTION = 'A flag to indicate whether to use XGBoost as the model for mokapots semi-supervised training loop. If set to False, mokapot will use its Percolator-like model.'
TEST_FDR_DESCRIPTION = """The false-discovery rate threshold at which to evaluate
        the learned models."""
FOLDS_DESCRIPTION = """The number of cross-validation folds to use. PSMs originating
        from the same mass spectrum are always in the same fold."""
WORKERS_DESCRIPTION = """The number of processes to use for model training. More workers
        will require more memory, but will typically decrease the total
        run time. An integer exceeding the number of folds will have
        no additional effect. Note that logging messages will be garbled
        if more than one worker is enabled."""
MAX_ITER_DESCRIPTION = 'The number of iterations to preform during the semi-supervized training loop'
TIMSCORE_DESCRIPTION = 'A flag to denote the use of PaSER SQT files. Setting this to True will 1) use timscores as a feature in mokapot and 2) output a PaSER compatible DTASelect-filter.txt file.'
MAX_SEED_DESCRIPTION = 'The random seed to use for reproducibility.'

# IP2Mokapot arguments
MASS_ALIGNMENT_DESCRIPTION = 'Alignment of masses within each SQT file, with recalculation of ppm based on a N-dimensional function fit to the mass-ppm drift versus retention time relationship (where N is specified by --mass_alignment_dim)'
MASS_ALIGNMENT_DIM_DESCRIPTION = 'The number of dimensions to use in the alignment function, with 1 being linear and 2 being quadratic, etc.'
MASS_ALIGNMENT_PERCENTILE_DESCRIPTION = 'The xcorr percentile for which to fit the mass alignment function, represented as an integer.'

MAX_MLINE_DESCRIPTION = 'The maximum number of m lines to be used as input for mokapot per spectra.'
XCORR_FILTER_DESCRIPTION = 'The xcorr value to be used as a filter for PSMs before they are included in the semi-supervised training loop, only PSMs with xcorr values greater than this will be included.'

VERBOSITY_DESCRIPTION = 'The level of verbosity to be used, with options ranging from 0 (ERROR), 1 (WARNING), 2 (INFO), to 3 (DEBUG).'
FILTER_LEVEL_DESCRIPTION = 'Corresponds to the -t value in the dtaselect.params file and determines the level of filtering for duplicate PSMs, with options ranging from 1 (keep all), 2 (keep best-scoring PSMs after grouping by sequence and charge), to 3 (keep best-scoring PSMs after grouping by only sequence).'

# Help msg
HELP_DESCRIPTION = """

Welcome to MokaFilter! This tutorial will help you get started with our app.

Start by uploading one or multiple SQT files. These files contain the raw peptide search matches (PSMs) from the Prolucid search engine.
Then, upload one or more FASTA files. These files should be the ones used to generate the SQT files.
Upload the search.xml file to specify the Prolucid search parameters.
Lastly, upload the dtaselect.params file which holds the filtering parameters.
Once all necessary files are uploaded, click "run" and MokaFilter will generate a filtered results list of proteins, peptides, and spectra at the provided false discovery rate (FDR).

Note: MokaFilter uses the mokapot software:

Fondrie W. E. & Noble W. S. mokapot: Fast and Flexible Semisupervised Learning for Peptide Detection. J Proteome Res (2021) doi: 10.1021/acs.jproteome.0c01010. PMID: 33596079.
"""

EMPTY_DTA_FILTER = """DTASelect v2.1.12
/data/
/data/
ProLuCID 1.4 in SQT format.
 -p 1 -y 1 --trypstat --pfp 0.01 --modstat --extra --pI --DB --dm -t 1 --after HYMFWL --brief --quiet
true    Use criteria
0.0     Minimum peptide probability
0.01    Peptide global false discovery rate
0.0     Minimum protein probability
1.0     Protein false discovery rate
1       Minimum charge state
50      Maximum charge state
-1.0    Minimum ion proportion
10000   Maximum Sp rank
-1.0    Minimum Sp score
true    Define delta mass with respect to nearest isotope
Include Modified peptide inclusion
Half    Tryptic status requirement
false   Multiple, ambiguous IDs allowed
Ignore  Peptide validation handling
Salt step       Purge duplicate peptides by protein
false   Include only loci with unique peptide
true    Remove subset proteins
Ignore  Locus validation handling
0       Minimum modified peptides per locus
1       Minimum peptides per locus

Locus   Sequence Count  Spectrum Count  Sequence Coverage       Length  MolWt   pI      Validation Status       NSAF    EMPAI   Descriptive Name
Unique  FileName        XCorr   DeltCN  Conf%   M+H+    CalcM+H+        PPM     TotalIntensity  SpR     Prob Score      pI      IonProportion   Redundancy      Sequence        PTMIndex        PTMIndex Protein List
        Proteins        Peptide IDs     Spectra
Unfiltered      0      0      0
Filtered        0       0       0
Forward matches 0       0       0
Redundant Forward matches       0       0       0
Decoy matches   0       0       0
Redundant Decoy matches 0       0       0
Forward FDR     0.0     0.0     0.0

Classification  Nonredundant Proteins   Redundant Proteins
Unclassified    0       0  
"""