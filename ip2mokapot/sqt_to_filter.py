from __future__ import annotations

import argparse
from io import StringIO, TextIOWrapper
from pathlib import Path
from collections import Counter

import mokapot
import pandas as pd
from serenipy import dtaselectfilter
from xgboost import XGBClassifier

from .parsing import convert_to_csv, convert_to_moka, get_filter_results_moka
from .util import xml_to_dict, read_fasta, _parse_fasta_files, _parse_protein


# TODO: Add routine for running multiple sqt files independently/combined/grouped/ungrouped
# TODO: Fix wierd TextIO vs TextIOWrapper typing issue
# TODO: Make the ct-erm/n-term enzyme arg work (currently it only supports c-term)
# TODO: Make option to save intermediate mokapot files somewhere

def parse_args() -> argparse.Namespace:
    _parser = argparse.ArgumentParser(description='Arguments for Percolator to DtaSelect-Filter')
    _parser.add_argument('--sqts', required=True, action='append', type=str, help='path to SQT file')
    _parser.add_argument('--fastas', required=True, action='append', type=str,
                         help='path to FASTA file (must contain decoys and target proteins)')
    _parser.add_argument('--out', required=True, type=str, help='path to write output file')

    _parser.add_argument('--protein_fdr', required=False, type=float, default=0.01, help='protein level FDR')
    _parser.add_argument('--peptide_fdr', required=False, type=float, default=0.01, help='peptide level FDR')
    _parser.add_argument('--psm_fdr', required=False, type=float, default=0.01, help='psm level FDR')
    _parser.add_argument('--min_peptides', required=False, default=1, type=int,
                         help='min number of peptides required for a protein to be identified')

    _parser.add_argument('--search_xml', required=False,
                         help='path to the search.xml file. USing this will override the following variables: '
                              'enzyme_regex, missed_cleavage, min_length, semi')
    _parser.add_argument('--enzyme_regex', required=False, default='[KR]',
                         help='The regex string to determine enzyme sites')
    _parser.add_argument('--missed_cleavage', required=False, default=0, type=int,
                         help='Number of internal peptide missed cleavages')
    _parser.add_argument('--min_length', required=False, default=6, type=float, help='min peptide length')
    _parser.add_argument('--max_length', required=False, default=50, type=float, help='max peptide length')
    _parser.add_argument('--semi', required=False, default=False, type=bool,
                         help='flag for whether peptides are semi enzymatic. set to True if both ends follow adhere '
                              'to the enzyme_regex')
    _parser.add_argument('--decoy_prefix', required=False, default='Reverse_', type=str,
                         help='decoy prefix found int FASTA and SQT files. For IP2 use Reverse_')

    _parser.add_argument('--xgboost', required=False, default=True, type=bool,
                         help='Use Xbgoost. If false will use percolator like algorithm')
    _parser.add_argument('--test_fdr', required=False, default=0.01, type=float,
                         help='test FDR to use during semi-supervized training loop')
    _parser.add_argument('--folds', required=False, default=3, type=int, help='number of K-Folds to preform')
    _parser.add_argument('--workers', required=False, default=1, type=int,
                         help='number of workers (threads) to use for semi-supervized training loop')
    _parser.add_argument('--max_iter', required=False, default=10, type=int,
                         help='number of iterations to preform during the semi-supervized training loop')

    return _parser.parse_args()


def run():
    """
    mokafilter entrypoint in setup.py.
    Parse args and convert sqts and fasta files into TextIO
    """
    args = parse_args()
    sqts = [Path(sqt).open() for sqt in args.sqts]
    sqt_stems = [str(Path(sqt).stem) for sqt in args.sqts]
    fastas = [Path(fasta).open() for fasta in args.fastas]
    search_xml = Path(args.search_xml).open()

    dta_filter_content = mokafilter(sqts, fastas, args.protein_fdr, args.peptide_fdr, args.psm_fdr, args.min_peptides,
                                    search_xml, args.enzyme_regex, args.missed_cleavage,
                                    args.min_length, args.max_length, args.semi, args.decoy_prefix, args.xgboost,
                                    args.test_fdr, args.folds, args.workers, sqt_stems, args.max_iter)

    with open(Path(args.out), 'w') as file:
        file.write(dta_filter_content)


def mokafilter(sqts: list[TextIOWrapper | StringIO], fastas: list[TextIOWrapper | StringIO], protein_fdr: float,
               peptide_fdr: float, psm_fdr: float, min_peptides: int, search_xml: TextIOWrapper | StringIO,
               enzyme_regex: str, missed_cleavage: int, min_length: int, max_length: int, semi: bool,
               decoy_prefix: str, xgboost: bool, test_fdr: float, folds: int, workers: int, sqt_stems: list[str],
               max_iter: int) -> str:
    """
    What a mess of code...

    Entrypoint for both CLI tool and streamlit app, as such all files but be of IO type (StringIO or TextIO)
    :return: str - the string contents of the output DTASelect-filter.txt file
    """

    sqt_dfs = [convert_to_csv(sqt, sqt_stem) for sqt, sqt_stem in zip(sqts, sqt_stems)]
    sqt_df = pd.concat(sqt_dfs, ignore_index=True)
    sqt_df = sqt_df[sqt_df['xcorr'] > 0.0]
    pin_df = convert_to_moka(sqt_df)

    # Override default/inputted attributes if search.xml file is not None
    if search_xml:
        xml_dict = xml_to_dict(search_xml)
        missed_cleavage = int(xml_dict['enzyme_info']['max_num_internal_mis_cleavage'])
        semi = int(xml_dict['enzyme_info']['specificity']) != 2
        enzyme_regex = f"[{''.join(xml_dict['enzyme_info']['residues']['residue'])}]"
        enzyme_term = 'cterm' if xml_dict['enzyme_info']['type'] == 'true' else 'nterm'
        min_length = int(xml_dict['peptide_length_limits']['minimum'])

    fasta_elems = [_parse_protein(entry) for entry in _parse_fasta_files(fastas)]
    fasta_dict = {e[0]: {'sequence': e[1], 'description': e[2]} for e in fasta_elems}

    psms = mokapot.read_pin(pin_files=pin_df)

    # Slightly modified version of read_fasta which doesn't require opening files
    proteins = read_fasta(fasta=fasta_elems,
                          enzyme=enzyme_regex,
                          missed_cleavages=missed_cleavage,
                          min_length=min_length,
                          max_length=max_length,
                          semi=semi,
                          decoy_prefix=decoy_prefix)
    psms.add_proteins(proteins)

    if xgboost:
        estimator = mokapot.Model(XGBClassifier(objective='binary:logistic', nthread=4, seed=42), max_iter=max_iter)
    else:
        estimator = mokapot.PercolatorModel(max_iter=max_iter)

    results, models = mokapot.brew(psms=psms,
                                   model=estimator,
                                   test_fdr=test_fdr,
                                   folds=folds,
                                   max_workers=workers)

    # results.to_txt(dest_dir=str(sqt_path.parent), file_root=str(sqt_path.stem))
    print(models)

    # separate protein, peptide, and psm dataframes
    target_psm_results, target_peptide_results, target_protein_results = results.confidence_estimates['psms'], \
        results.confidence_estimates['peptides'], results.confidence_estimates['proteins']
    decoy_psm_results, decoy_peptide_results, decoy_protein_results = results.decoy_confidence_estimates['psms'], \
        results.decoy_confidence_estimates['peptides'], results.decoy_confidence_estimates['proteins']

    # Filter each dataframe according the specified level FDR value
    filtered_target_psm_results = target_psm_results[target_psm_results['mokapot q-value'] <= psm_fdr]
    filtered_target_peptide_results = target_peptide_results[target_peptide_results['mokapot q-value'] <= peptide_fdr]
    filtered_target_protein_results = target_protein_results[target_protein_results['mokapot q-value'] <= protein_fdr]
    filtered_decoy_psm_results = decoy_psm_results[decoy_psm_results['mokapot q-value'] <= psm_fdr]
    filtered_decoy_peptide_results = decoy_peptide_results[decoy_peptide_results['mokapot q-value'] <= peptide_fdr]
    filtered_decoy_protein_results = decoy_protein_results[decoy_protein_results['mokapot q-value'] <= protein_fdr]

    # Convert mokapot result dataframes into serenipy DTASelect-filter results
    # Keep only the results which have >= min peptide lines
    target_filter_results = get_filter_results_moka(sqt_df, filtered_target_psm_results,
                                                    filtered_target_peptide_results,
                                                    filtered_target_protein_results, fasta_dict)
    target_filter_results = [result for result in target_filter_results if len(result.peptide_lines) >= min_peptides]
    decoy_filter_results = get_filter_results_moka(sqt_df, filtered_decoy_psm_results, filtered_decoy_peptide_results,
                                                   filtered_decoy_protein_results, fasta_dict)
    decoy_filter_results = [result for result in decoy_filter_results if len(result.peptide_lines) >= min_peptides]
    filter_results = target_filter_results + decoy_filter_results

    # Assign the unique property to peptide lines
    peptide_sets = [{line.sequence for line in result.peptide_lines} for result in filter_results]
    peptide_counts = Counter([peptide for peptide_set in peptide_sets for peptide in peptide_set])
    for result in filter_results:
        for peptide_line in result.peptide_lines:
            peptide_line.unique = '*' if peptide_counts[peptide_line.sequence] == 1 else ''

    # Merge DTASelect-filter results for proteins which mokapot has grouped together. proteins seem to be grouped only
    # if the share the same peptide lines
    peptide_groups = [g.split(', ') for g in filtered_target_protein_results['mokapot protein group'].values if
                      len(g.split(', ')) > 1]
    for peptide_group in peptide_groups:
        results = []
        for i, result in reversed(list(enumerate(filter_results))):
            if result.protein_lines[0].locus_name in peptide_group:
                results.append(result)
                filter_results.pop(i)
        if results:
            new_result = dtaselectfilter.DTAFilterResult(protein_lines=[result.protein_lines[0] for result in results],
                                                         peptide_lines=results[0].peptide_lines)
            filter_results.append(new_result)

    # Finally sort all filter results by sequence coverage
    filter_results.sort(key=lambda x: x.protein_lines[0].sequence_coverage, reverse=True)

    # Add protein line NSAF:
    safs = [protein_line.spectrum_count / protein_line.length for result in filter_results for protein_line in
            result.protein_lines]
    sum_safs, i = sum(safs), 0
    for result in filter_results:
        for protein_line in result.protein_lines:
            protein_line.nsaf = safs[i] / sum_safs
            i += 1

    # Add peptide line redundancy
    peptide_counts = Counter(
        [peptide_line.sequence for result in filter_results for peptide_line in result.peptide_lines])
    for result in filter_results:
        for peptide_line in result.peptide_lines:
            peptide_line.redundancy = peptide_counts[peptide_line.sequence]

    # Finalize DTASelect-filter.txt lines
    h_lines = ['DTASelect v2.1.12\n',
               '\n',
               'Locus	Sequence Count	Spectrum Count	Sequence Coverage	Length	MolWt	pI	Validation Status	NSAF	EMPAI	Descriptive Name	HRedundancy	LRedundancy	MRedundancy\n',
               'Unique	FileName	XCorr	DeltCN	Conf%	M+H+	CalcM+H+	PPM	TotalIntensity	SpR	Prob Score	pI	IonProportion	Redundancy	Sequence	RetTime	PTMIndex	PTMIndex Protein List\n'
               ]

    unfiltered_proteins = len(target_protein_results) + len(decoy_protein_results)
    unfiltered_peptides = len(target_peptide_results) + len(decoy_peptide_results)
    unfiltered_psms = len(target_psm_results) + len(decoy_psm_results)

    filtered_proteins = len(filtered_target_protein_results) + len(filtered_decoy_protein_results)
    filtered_peptides = len(filtered_target_peptide_results) + len(filtered_decoy_peptide_results)
    filtered_psms = len(filtered_target_psm_results) + len(filtered_decoy_psm_results)

    forward_proteins = len(filtered_target_protein_results)
    forward_peptides = len(filtered_target_peptide_results)
    forward_psms = len(filtered_target_psm_results)

    decoy_proteins = len(filtered_decoy_protein_results)
    decoy_peptides = len(filtered_decoy_peptide_results)
    decoy_psms = len(filtered_decoy_psm_results)

    protein_fdr = len(filtered_decoy_protein_results) / (len(filtered_target_protein_results) +
                                                         len(filtered_decoy_protein_results)) * 100
    peptide_fdr = len(filtered_decoy_peptide_results) / (len(filtered_target_peptide_results) +
                                                         len(filtered_decoy_peptide_results)) * 100
    psm_fdr = len(filtered_decoy_psm_results) / (len(filtered_target_psm_results) +
                                                 len(filtered_decoy_psm_results)) * 100

    end_lines = [f'	Proteins	Peptide IDs	Spectra\n',
                 f'Unfiltered	{unfiltered_proteins}  {unfiltered_peptides}  {unfiltered_psms}\n',
                 f'Filtered	{filtered_proteins}  {filtered_peptides}  {filtered_psms}\n',
                 f'Forward	{forward_proteins}  {forward_peptides}  {forward_psms}\n',
                 f'Redundant Forward matches	{0}  {0}  {0}\n',
                 f'Decoy matches	{decoy_proteins}  {decoy_peptides}  {decoy_psms}\n',
                 f'Redundant Decoy matches	{0}  {0}  {0}\n',
                 f'Forward FDR	{protein_fdr}	{peptide_fdr}	{psm_fdr}\n']

    dta_filter_content = dtaselectfilter.to_dta_select_filter(
        version=dtaselectfilter.DtaSelectFilterVersion.V2_1_12_paser,
        h_lines=h_lines,
        dta_filter_results=filter_results,
        end_lines=end_lines)

    return dta_filter_content


if __name__ == '__main__':
    run()
