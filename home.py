from io import StringIO
import streamlit as st
from ip2mokapot.sqt_to_filter import mokafilter

st.title('IP2MokaFilter')
st.text('Converts SQT to DTASelect-filter.txt with mokapot')


sqts = st.file_uploader(label='sqt', type='.sqt', accept_multiple_files=True)
fastas = st.file_uploader(label='fasta', type='.fasta', accept_multiple_files=True)
search_xml = st.file_uploader(label='search_xml', type='.xml')

with st.expander('Advanced'):
    protein_fdr = st.number_input(label='protein_fdr', value=0.01)
    peptide_fdr = st.number_input(label='peptide_fdr', value=0.01)
    psm_fdr = st.number_input(label='psm_fdr', value=0.01)
    min_peptides = st.number_input(label='min_peptides', value=1)

    enzyme_regex = st.text_input(label='enzyme_regex', value='[KR]')
    enzyme_term = st.checkbox(label='enzyme_term', value=True)

    missed_cleavage = st.number_input(label='missed_cleavage', value=1)
    min_length = st.number_input(label='min_length', value=6)
    max_length = st.number_input(label='max_length', value=50)
    semi = st.checkbox(label='semi', value=False)
    decoy_prefix = st.text_input(label='decoy_prefix', value='Reverse_')
    xgboost = st.checkbox(label='xgboost', value=True)
    test_fdr = st.number_input(label='test_fdr', value=0.01)
    folds = st.number_input(label='folds', value=3)
    workers = st.number_input(label='workers', value=1)
    max_itr = st.number_input(label='max_itr', value=10)

if st.button('start') and sqts and fastas:

    sqt_ios = [StringIO(sqt.getvalue().decode("utf-8")) for sqt in sqts]
    fasta_ios = [StringIO(fasta.getvalue().decode("utf-8")) for fasta in fastas]
    sqt_stems = [''.join(sqt.name.split('.')[:-1]) for sqt in sqts]

    search_xml_io = None
    if search_xml:
        search_xml_io = StringIO(search_xml.getvalue().decode("utf-8"))

    dta_filter_content = mokafilter(sqt_ios, fasta_ios, protein_fdr, peptide_fdr, psm_fdr, min_peptides,
                                    search_xml_io, enzyme_regex, enzyme_term, missed_cleavage, min_length, max_length, semi,
                                    decoy_prefix, xgboost, test_fdr, folds, workers, sqt_stems, max_itr)

    st.download_button(label='Download DTASelect-filter.txt', data=dta_filter_content, file_name='DTASelect-filter.txt')
