from io import StringIO
import streamlit as st
from ip2mokapot.filter import mokafilter
from ip2mokapot.config import *

st.title('MokaFilter 	:coffee:')
st.text('Generates DTASelect-filter.txt with mokapot')

with st.expander('Help'):
    st.markdown(HELP_DESCRIPTION)

sqts = st.file_uploader(label='sqt', type='.sqt', accept_multiple_files=True, help=SQTS_DESCRIPTION)
fastas = st.file_uploader(label='fasta', type='.fasta', accept_multiple_files=True, help=FASTAS_DESCRIPTION)
search_xml = st.file_uploader(label='search_xml', type='.xml', help=SEARCH_XML_DESCRIPTION)

c1, c2 = st.columns(2)
dta_params_input_type = c1.radio(
    "dta_params_input",
    ('text','file'))
dta_params, dta_params_txt = None, None
if dta_params_input_type == 'file':
    dta_params = c2.file_uploader(label='dta_params', type='.params', help=DTASELECT_PARAMS_DESCRIPTION)
else:
    dta_params_txt = c2.text_input(label='dta_params', value='', help=DTASELECT_PARAMS_DESCRIPTION)


with st.expander('Advanced'):
    protein_fdr, peptide_fdr, psm_fdr, min_peptides = 0.01, 0.01, 0.01, 1
    if not dta_params:
        protein_fdr = st.number_input(label='protein_fdr', value=0.01, help=PROTEIN_FDR_DESCRIPTION)
        peptide_fdr = st.number_input(label='peptide_fdr', value=0.01, help=PEPTIDE_FDR_DESCRIPTION)
        psm_fdr = st.number_input(label='psm_fdr', value=0.01, help=PSM_FDR_DESCRIPTION)
        min_peptides = st.number_input(label='min_peptides', value=1, help=MIN_PEPTIDES_DESCRIPTION)

    missed_cleavage, semi, enzyme_regex, enzyme_term, min_length = 1, False, '[KR]', True, 6
    if not search_xml:
        missed_cleavage = st.number_input(label='missed_cleavage', value=1, help=MISSED_CLEAVAGE_DESCRIPTION)
        semi = st.checkbox(label='semi', value=False, help=SEMI_DESCRIPTION)
        enzyme_regex = st.text_input(label='enzyme_regex', value='[KR]', help=ENZYME_REGEX_DESCRIPTION)
        enzyme_term = st.checkbox(label='enzyme_term', value=True, help=ENZYME_TERM_DESCRIPTION)
        min_length = st.number_input(label='min_length', value=6, help=MIN_LENGTH_DESCRIPTION)

    max_length = st.number_input(label='max_length', value=50, help=MAX_LENGTH_DESCRIPTION)
    decoy_prefix = st.text_input(label='decoy_prefix', value='Reverse_', help=DECOY_PREFIX_DESCRIPTION)
    xgboost = st.checkbox(label='xgboost', value=True, help=XGBOOST_DESCRIPTION)
    test_fdr = st.number_input(label='test_fdr', value=0.01, help=TEST_FDR_DESCRIPTION)
    folds = st.number_input(label='folds', value=3, help=FOLDS_DESCRIPTION)
    workers = st.number_input(label='workers', value=1, help=WORKERS_DESCRIPTION)
    max_itr = st.number_input(label='max_itr', value=10, help=MAX_ITER_DESCRIPTION)

    timscore = st.checkbox(label='timscore', value=False, help=TIMSCORE_DESCRIPTION)
    mass_alignment = st.checkbox(label='mass_alignment', value=True, help=MASS_ALIGNMENT_DESCRIPTION)
    mass_alignment_dim = st.number_input(label='mass_alignment_dim', value=1, help=MASS_ALIGNMENT_DIM_DESCRIPTION)
    mass_alignment_percentile = st.number_input(label='mass_alignment_percentile', value=95, help=MASS_ALIGNMENT_PERCENTILE_DESCRIPTION)

    max_mline = st.number_input(label='max_mline', value=5, help=MAX_MLINE_DESCRIPTION)
    xcorr_filter = st.number_input(label='xcorr_filter', value=0.0, help=XCORR_FILTER_DESCRIPTION)

    use_random_seed = st.checkbox(label='random_seed', value=False, help=MAX_SEED_DESCRIPTION)

    random_seed = None
    if use_random_seed is True:
        random_seed = st.number_input(label='random_seed', value=42)

if st.button('Run'):

    if not sqts or not fastas:
        st.warning('Please upload the required file types: SQT & FASTA')
        st.stop()

    sqt_ios = [StringIO(sqt.getvalue().decode("utf-8")) for sqt in sqts]
    fasta_ios = [StringIO(fasta.getvalue().decode("utf-8")) for fasta in fastas]
    sqt_stems = [''.join(sqt.name.split('.')[:-1]) for sqt in sqts]

    search_xml_io = None
    if search_xml:
        search_xml_io = StringIO(search_xml.getvalue().decode("utf-8"))

    dta_params_io = None
    if dta_params:
        dta_params_io = StringIO(dta_params.getvalue().decode("utf-8"))
    elif dta_params_txt:
        dta_params_io = StringIO(dta_params_txt)

    alignment_figs, dta_filter_content = mokafilter(sqt_ios, fasta_ios, protein_fdr, peptide_fdr, psm_fdr, min_peptides,
                                    search_xml_io, enzyme_regex, enzyme_term, missed_cleavage, min_length, max_length,
                                    semi,
                                    decoy_prefix, xgboost, test_fdr, folds, workers, sqt_stems, max_itr, timscore,
                                    mass_alignment,
                                    max_mline, random_seed, dta_params_io, xcorr_filter, mass_alignment_dim, mass_alignment_percentile)

    if alignment_figs:
        for fig, stem in zip(alignment_figs, sqt_stems):
            st.caption(stem)
            st.pyplot(fig)

    st.download_button(label='Download DTASelect-filter.txt', data=dta_filter_content.read(), file_name='DTASelect-filter.txt')
