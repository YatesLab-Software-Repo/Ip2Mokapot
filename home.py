from io import StringIO
import streamlit as st
from ip2mokapot.sqt_to_filter import mokafilter

st.title('IP2MokaFilter')
st.text('Converts SQT to DTASelect-filter.txt with mokapot')


sqts = st.file_uploader(label='sqt', type='.sqt', accept_multiple_files=True)
fastas = st.file_uploader(label='fasta', type='.fasta', accept_multiple_files=True)
search_xml = st.file_uploader(label='search_xml', type='.xml')
dta_params = st.file_uploader(label='dta_params', type='.params')

with st.expander('Advanced'):
    protein_fdr = st.number_input(label='protein_fdr', value=0.01)
    peptide_fdr = st.number_input(label='peptide_fdr', value=0.01)
    psm_fdr = st.number_input(label='psm_fdr', value=0.01)
    min_peptides = st.number_input(label='min_peptides', value=1)

    missed_cleavage, semi, enzyme_regex, enzyme_term, min_length = None, None, None, None, None
    if not search_xml:
        missed_cleavage = st.number_input(label='missed_cleavage', value=1)
        semi = st.checkbox(label='semi', value=False)
        enzyme_regex = st.text_input(label='enzyme_regex', value='[KR]')
        enzyme_term = st.checkbox(label='enzyme_term', value=True)
        min_length = st.number_input(label='min_length', value=6)

    max_length = st.number_input(label='max_length', value=50)
    decoy_prefix = st.text_input(label='decoy_prefix', value='Reverse_')
    xgboost = st.checkbox(label='xgboost', value=True)
    test_fdr = st.number_input(label='test_fdr', value=0.01)
    folds = st.number_input(label='folds', value=3)
    workers = st.number_input(label='workers', value=1)
    max_itr = st.number_input(label='max_itr', value=10)

    timscore = st.checkbox(label='timscore', value=False)
    mass_alignment = st.checkbox(label='mass_alignment', value=True)
    max_mline = st.number_input(label='max_mline', value=5)
    use_random_seed = st.checkbox(label='random_seed', value=False)

    random_seed= None
    if use_random_seed is True:
        random_seed = st.number_input(label='random_seed', value=42)

if st.button('start') and sqts and fastas:

    sqt_ios = [StringIO(sqt.getvalue().decode("utf-8")) for sqt in sqts]
    fasta_ios = [StringIO(fasta.getvalue().decode("utf-8")) for fasta in fastas]
    sqt_stems = [''.join(sqt.name.split('.')[:-1]) for sqt in sqts]

    search_xml_io = None
    if search_xml:
        search_xml_io = StringIO(search_xml.getvalue().decode("utf-8"))

    dta_params_io = None
    if dta_params:
        dta_params_io = StringIO(dta_params.getvalue().decode("utf-8"))


    dta_filter_content = mokafilter(sqt_ios, fasta_ios, protein_fdr, peptide_fdr, psm_fdr, min_peptides,
                                    search_xml_io, enzyme_regex, enzyme_term, missed_cleavage, min_length, max_length, semi,
                                    decoy_prefix, xgboost, test_fdr, folds, workers, sqt_stems, max_itr, timscore, mass_alignment,
                                    max_mline, random_seed, dta_params_io)

    st.download_button(label='Download DTASelect-filter.txt', data=dta_filter_content, file_name='DTASelect-filter.txt')
