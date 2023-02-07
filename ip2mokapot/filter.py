import argparse
import logging
from io import StringIO
from pathlib import Path
from collections import Counter
from typing import List, Union, Tuple, Iterator

import matplotlib
import numpy as np
import pandas as pd
import mokapot
from mokapot.parsers.fasta import _parse_protein, _parse_fasta_files, read_fasta
from tabulate import tabulate
from serenipy import dtaselectfilter
from typing.io import IO
from xgboost import XGBClassifier

from .parsing import convert_to_csv, convert_to_moka, get_filter_results_moka, align_mass, parse_dta_args
from .util import xml_to_dict
from .config import *


# TODO: Make option to save intermediate mokapot files somewhere
# TODO: Add option to train hyper params
# TODO: Make Protein digestion threaded (since its slow for semi)
# TODO: Fix enzyme regex and tryptic status issue, For complex regexes the tryptic status will not work right

def float_range(mini,maxi):
    """Return function handle of an argument type function for
       ArgumentParser checking a float range: mini <= arg <= maxi
         mini - minimum acceptable argument
         maxi - maximum acceptable argument"""

    # Define the function with default arguments
    def float_range_checker(arg):
        """New Type function for argparse - a float within predefined range."""

        try:
            f = float(arg)
        except ValueError:
            raise argparse.ArgumentTypeError("must be a floating point number")
        if f < mini or f > maxi:
            raise argparse.ArgumentTypeError("must be in range [" + str(mini) + " .. " + str(maxi)+"]")
        return f

    # Return function handle to checking function
    return float_range_checker
def parse_args() -> argparse.Namespace:
    """
    Argument parser for ip2mokapot
    :return: argparse.Namespace - parsed arsg
    """
    _parser = argparse.ArgumentParser(description='Arguments for MokaFilter',
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    _parser.add_argument('--sqts', required=True, nargs='+', type=str, help=SQTS_DESCRIPTION)
    _parser.add_argument('--fastas', required=True, nargs='+', type=str, help=FASTAS_DESCRIPTION)
    _parser.add_argument('--out', required=True, type=str, help=OUT_DESCRIPTION)
    _parser.add_argument('--search_xml', required=False, type=str, default=None, help=SEARCH_XML_DESCRIPTION)
    _parser.add_argument('--dta_params', required=False, type=str, default=None, help=DTASELECT_PARAMS_DESCRIPTION)

    _parser.add_argument('--protein_fdr', required=False, type=float, default=0.01, help=PROTEIN_FDR_DESCRIPTION)
    _parser.add_argument('--peptide_fdr', required=False, type=float, default=0.01, help=PEPTIDE_FDR_DESCRIPTION)
    _parser.add_argument('--psm_fdr', required=False, type=float, default=0.01, help=PSM_FDR_DESCRIPTION)
    _parser.add_argument('--min_peptides', required=False, default=1, type=int, help=MIN_PEPTIDES_DESCRIPTION)

    _parser.add_argument('--enzyme_regex', required=False, default='[KR]', help=ENZYME_REGEX_DESCRIPTION)
    _parser.add_argument('--enzyme_term', required=False, default=True, help=ENZYME_TERM_DESCRIPTION)
    _parser.add_argument('--missed_cleavage', required=False, default=0, type=int, help=MISSED_CLEAVAGE_DESCRIPTION)
    _parser.add_argument('--min_length', required=False, default=6, type=float, help=MIN_LENGTH_DESCRIPTION)
    _parser.add_argument('--max_length', required=False, default=50, type=float, help=MAX_LENGTH_DESCRIPTION)
    _parser.add_argument('--semi', required=False, default=False, action=argparse.BooleanOptionalAction,
                         help=SEMI_DESCRIPTION)
    _parser.add_argument('--decoy_prefix', required=False, default='Reverse_', type=str, help=DECOY_PREFIX_DESCRIPTION)

    _parser.add_argument('--xgboost', required=False, default=False, action=argparse.BooleanOptionalAction,
                         help=XGBOOST_DESCRIPTION)
    _parser.add_argument('--test_fdr', required=False, default=0.01, type=float, help=TEST_FDR_DESCRIPTION)
    _parser.add_argument('--folds', required=False, default=3, type=int, help=FOLDS_DESCRIPTION)
    _parser.add_argument('--workers', required=False, default=1, type=int, help=WORKERS_DESCRIPTION)
    _parser.add_argument('--max_iter', required=False, default=10, type=int, help=MAX_ITER_DESCRIPTION)

    _parser.add_argument('--timscore', default=False, action=argparse.BooleanOptionalAction, help=TIMSCORE_DESCRIPTION)
    _parser.add_argument('--mass_alignment', default=True, action=argparse.BooleanOptionalAction,
                         help=MASS_ALIGNMENT_DESCRIPTION)
    _parser.add_argument('--mass_alignment_dim', required=False, default=1, type=int,
                         help=MASS_ALIGNMENT_DIM_DESCRIPTION)
    _parser.add_argument('--mass_alignment_percentile', required=False, default=95, type=int,
                         help=MASS_ALIGNMENT_PERCENTILE_DESCRIPTION)

    _parser.add_argument('--max_mline', required=False, default=None, type=int, help=MAX_MLINE_DESCRIPTION)
    _parser.add_argument('--seed', required=False, default=None, type=int, help=MAX_SEED_DESCRIPTION)
    _parser.add_argument('--xcorr_filter', required=False, default=None, type=float, help=XCORR_FILTER_DESCRIPTION)
    _parser.add_argument('--verbosity', required=False, default=1, type=int, choices=[0,1,2,3], help=VERBOSITY_DESCRIPTION)
    _parser.add_argument('--filter_level', required=False, default=0, type=int, choices=[0,1,2], help=FILTER_LEVEL_DESCRIPTION)

    return _parser.parse_args()


def io_gen(paths: List[Path], ext) -> Iterator[IO[str]]:
    for path in paths:
        if path.is_file():
            with path.open() as file:
                yield file
        else:
            if path.is_dir():
                pathlist = path.rglob(f'*{ext}')
                for p in pathlist:
                    with p.open() as file:
                        yield file


def stem_gen(paths: List[Path], ext) -> Iterator[str]:
    for path in paths:
        if path.is_file():
            yield path.stem
        else:
            if path.is_dir():
                pathlist = path.rglob(f'*{ext}')
                for p in pathlist:
                    yield p.stem


def run():
    """
    mokafilter entrypoint in setup.py.
    Parse args and convert sqts and fasta files into TextIO

    """

    args = parse_args()

    verbosity_dict = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }

    logging.basicConfig(
        format=("[{levelname}] {message}"),
        style="{",
        level=verbosity_dict[args.verbosity],
    )

    sqt_paths = [Path(sqt) for sqt in args.sqts]
    sqt_ios = io_gen(sqt_paths, '.sqt')
    sqt_stems = list(stem_gen(sqt_paths, '.sqt'))

    fasta_paths = [Path(fasta) for fasta in args.fastas]
    fasta_ios = io_gen(fasta_paths, '.fasta')
    fasta_stems = list(stem_gen(fasta_paths, '.fasta'))

    search_xml = None
    if args.search_xml:
        search_xml = Path(args.search_xml).open()

    dta_params = None
    if args.dta_params:
        dta_params = Path(args.dta_params).open()

    print(args)
    alignment_figs, pin_df, dta_filter_content = mokafilter((sqt_ios, sqt_stems), (fasta_ios, fasta_stems),
                                                            args.protein_fdr, args.peptide_fdr, args.psm_fdr,
                                                            args.min_peptides,
                                                            search_xml, args.enzyme_regex, args.enzyme_term,
                                                            args.missed_cleavage,
                                                            args.min_length, args.max_length, args.semi,
                                                            args.decoy_prefix, args.xgboost,
                                                            args.test_fdr, args.folds, args.workers, args.max_iter,
                                                            args.timscore,
                                                            args.mass_alignment, args.max_mline, args.seed, dta_params,
                                                            args.xcorr_filter,
                                                            args.mass_alignment_dim, args.mass_alignment_percentile,
                                                            args.filter_level)

    with open(Path(args.out), 'w') as file:
        file.write(dta_filter_content.read())

    pin_df.to_csv(args.out + '.pin', index=False)

    if alignment_figs:
        for fig, stem in zip(alignment_figs, sqt_stems):
            if fig is None:
                continue
            png_name = f'{stem}_mass_alignment.png'
            print(f'Saving alignment plot to {png_name}')
            fig.savefig(png_name)

    if dta_params:
        dta_params.close()

    if search_xml:
        search_xml.close()


def mokafilter(sqts: Tuple[Iterator[IO[str]], List[str]],
               fastas: Tuple[Iterator[IO[str]], List[str]],
               protein_fdr: float,
               peptide_fdr: float,
               psm_fdr: float,
               min_peptides: int,
               search_xml: Union[IO[str], None],
               enzyme_regex: str,
               enzyme_term: bool,
               missed_cleavage: int,
               min_length: int,
               max_length: int,
               semi: bool,
               decoy_prefix: str,
               xgboost: bool,
               test_fdr: float,
               folds: int, workers: int,
               max_iter: int,
               timscore: bool,
               mass_alignment: bool,
               max_mline: int,
               seed: Union[int, None],
               dta_params: Union[IO[str], None],
               xcorr_filter: float,
               mass_alignment_dim: int,
               mass_alignment_percentile: float,
               filter_level: int) -> (List[matplotlib.figure.Figure], pd.DataFrame, IO[str]):
    """
    What a mess of code...

    Entrypoint for both CLI tool and streamlit app, as such all files but be of IO type (StringIO or TextIO)
    :return: str - the string contents of the output DTASelect-filter.txt file
    """

    if dta_params:
        dta_args = parse_dta_args(dta_params.read().rstrip())
        fp_fdr = float(dta_args.get('--fp', 1.0))
        pfp_fdr = float(dta_args.get('--pfp', 1.0))
        sfp_fdr = float(dta_args.get('--sfp', 1.0))
        protein_fdr, peptide_fdr, psm_fdr = pfp_fdr, sfp_fdr, fp_fdr

        min_peptides = int(dta_args.get('-p', min_peptides))
        timscore = dta_args.get('--timscore', timscore)
        filter_level = int(dta_args.get('-t', filter_level))

    if search_xml:
        xml_dict = xml_to_dict(search_xml)
        missed_cleavage = int(xml_dict['enzyme_info']['max_num_internal_mis_cleavage'])
        if missed_cleavage == -1:
            missed_cleavage = max_length
        specificity = int(xml_dict['enzyme_info']['specificity'])
        semi = specificity != 2
        if specificity == 0:
            enzyme_regex = '.'
            missed_cleavage = max_length
        else:
            enzyme_regex = f"[{''.join(xml_dict['enzyme_info']['residues']['residue'])}]"
        enzyme_term = xml_dict['enzyme_info']['type'] == 'true'
        min_length = int(xml_dict['peptide_length_limits']['minimum'])

    tabulated_args = tabulate([
        ["sqt files", list(sqts[1])],
        ["fasta files", list(fastas[1])],
        ["search_xml", True if search_xml else False],
        ['dta_params', True if dta_params else False],
        ['protein_fdr', protein_fdr],
        ['peptide_fdr', peptide_fdr],
        ["psm_fdr", psm_fdr],
        ["min_peptides", min_peptides],
        ["enzyme_regex", enzyme_regex],
        ["enzyme_term", enzyme_term],
        ["missed_cleavage", missed_cleavage],
        ["min_length", min_length],
        ["max_length", max_length],
        ["semi", semi],
        ["decoy_prefix", decoy_prefix],
        ["xgboost", xgboost],
        ['test_fdr', test_fdr],
        ["folds", folds],
        ["workers", workers],
        ["max_iter", max_iter],
        ["timscore", timscore],
        ["mass_alignment", mass_alignment],
        ['mass_alignment_dim', mass_alignment_dim],
        ['mass_alignment_percentile', mass_alignment_percentile],
        ["max_mline", max_mline],
        ["seed", seed],
        ['xcorr_filter', xcorr_filter],
        ['filter_level', filter_level]
    ], headers=['Argument', 'Value'], missingval='None')

    print(tabulated_args)

    # Set the random seed:
    if seed:
        print(f'Setting random seed: {seed}')
        np.random.seed(seed)

    sqt_dfs = [convert_to_csv(sqt_io, sqt_stem, xcorr_filter, max_mline) for sqt_io, sqt_stem in zip(sqts[0], sqts[1])]
    alignment_figs = None
    if mass_alignment is True:
        print(f'Aligning masses...')
        figs = [align_mass(sqt_df, mass_alignment_dim, mass_alignment_percentile) for sqt_df in sqt_dfs]
        alignment_figs = figs

    sqt_df = pd.concat(sqt_dfs, ignore_index=True)
    pin_df = convert_to_moka(sqt_df)

    print(pin_df)

    fasta_elems = [_parse_protein(entry) for entry in _parse_fasta_files(fastas[0])]
    fasta_dict = {e[0]: {'sequence': e[1], 'description': e[2]} for e in fasta_elems}

    """if enzyme_term is True:
        pin_df['tryp'] = [
            int(peptide[0] in enzyme_regex[1:-1]) + int(strip_modifications(peptide[:-2])[-1] in enzyme_regex[1:-1]) for
            peptide in
            sqt_df['sequence']]
    else:
        pin_df['tryp'] = [int(peptide[2] in enzyme_regex[1:-1]) + int(peptide[-1] in enzyme_regex[1:-1]) for peptide in
                          sqt_df['sequence']]"""

    if any(sqt_df['tims_score']) and timscore is True:
        pin_df['tims_score'] = sqt_df['tims_score']
        pin_df['tims_score'] = pin_df['tims_score'].fillna(value=0)

    psms = mokapot.read_pin(pin_files=pin_df)

    # Slightly modified version of read_fasta which doesn't require opening files
    proteins = read_fasta(fasta=fasta_elems,
                          enzyme=enzyme_regex,
                          missed_cleavages=missed_cleavage,
                          min_length=min_length,
                          max_length=max_length,
                          semi=semi,
                          decoy_prefix=decoy_prefix,
                          enzyme_term=enzyme_term)

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
                                                    filtered_target_peptide_results, filtered_target_protein_results,
                                                    fasta_dict)

    decoy_filter_results = get_filter_results_moka(sqt_df, filtered_decoy_psm_results,
                                                   filtered_decoy_peptide_results, filtered_decoy_protein_results,
                                                   fasta_dict)

    filter_results = target_filter_results + decoy_filter_results

    filter_results = [result for result in filter_results if result.protein_lines[0].sequence_count >= min_peptides]

    for result in filter_results:
        result.filter(level=filter_level)

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
    unfiltered_proteins = len(target_protein_results) + len(decoy_protein_results)
    unfiltered_peptides = len(target_peptide_results) + len(decoy_peptide_results)
    unfiltered_psms = len(target_psm_results) + len(decoy_psm_results)

    target_results = [result for result in filter_results if
                      any(['Reverse_' not in protein_line.locus_name for protein_line in result.protein_lines])]
    target_protein_groups = len([result.protein_lines[0].locus_name for result in target_results])
    target_proteins = sum([len(result.protein_lines) for result in target_results])

    target_peptide_charge_pairs = [(peptide_line.charge, peptide_line.sequence[2:-2]) for result in target_results for
                                   peptide_line in
                                   result.peptide_lines]
    target_peptides = len(set(target_peptide_charge_pairs))
    total_target_peptides = len(target_peptide_charge_pairs)
    target_spectra = sum([result.protein_lines[0].spectrum_count for result in target_results])

    decoy_results = [result for result in filter_results if
                     all(['Reverse_' in protein_line.locus_name for protein_line in result.protein_lines])]
    decoy_protein_groups = len([result.protein_lines[0].locus_name for result in decoy_results])
    decoy_proteins = sum([len(result.protein_lines) for result in decoy_results])

    decoy_peptide_charge_pairs = [(peptide_line.charge, peptide_line.sequence[2:-2]) for result in decoy_results for
                                  peptide_line in
                                  result.peptide_lines]
    decoy_peptides = len(set(decoy_peptide_charge_pairs))
    total_decoy_peptides = len(decoy_peptide_charge_pairs)
    decoy_spectra = sum([result.protein_lines[0].spectrum_count for result in decoy_results])

    try:
        protein_fdr = round(decoy_protein_groups / (decoy_protein_groups + target_protein_groups) * 100, 4)
    except ZeroDivisionError:
        protein_fdr = 'NA'

    try:
        peptide_fdr = round(decoy_peptides / (decoy_peptides + target_peptides) * 100, 4)
    except ZeroDivisionError:
        peptide_fdr = 'NA'
    try:
        spectra_fdr = round(decoy_spectra / (decoy_spectra + target_spectra) * 100, 4)
    except ZeroDivisionError:
        spectra_fdr = 'NA'

    end_lines = [f'	Proteins	Peptide IDs	Spectra\n',
                 f'Unfiltered	{unfiltered_proteins}  {unfiltered_peptides}  {unfiltered_psms}\n',
                 f'Filtered	{decoy_protein_groups + target_protein_groups}  {target_peptides + decoy_peptides}  {target_spectra + decoy_spectra}\n',
                 f'Forward	{target_protein_groups}  {target_peptides}  {target_spectra}\n',
                 f'Redundant Forward matches	{target_proteins}  {total_target_peptides}  {target_spectra}\n',
                 f'Decoy matches	{decoy_protein_groups}  {decoy_peptides}  {decoy_spectra}\n',
                 f'Redundant Decoy matches	{decoy_proteins}  {total_decoy_peptides}  {decoy_spectra}\n',
                 f'Forward FDR	{protein_fdr}	{peptide_fdr}	{spectra_fdr}\n']

    if timscore is True:
        h_lines = ['DTASelect v2.1.13\n',
                   tabulated_args,
                   '\n',
                   'Locus	Sequence Count	Spectrum Count	Sequence Coverage	Length	MolWt	pI	'
                   'Validation Status	NSAF	EMPAI	Descriptive Name	HRedundancy	LRedundancy	MRedundancy\n',
                   'Unique	FileName	XCorr	DeltCN	Conf%	M+H+	CalcM+H+	PPM	TotalIntensity	SpR	'
                   'Prob Score	pI	IonProportion	Redundancy	Measured_IM_Value	Predicted_IM_Value	IM_Score	'
                   'Sequence	ExperimentalMz	Corrected1k0	IonMobility	RetTime	PTMIndex	PTMIndex Protein List\n'
                   ]
        dta_filter_content = dtaselectfilter.to_dta_select_filter(
            version=dtaselectfilter.DtaSelectFilterVersion.V2_1_13_timscore,
            h_lines=h_lines,
            dta_filter_results=filter_results,
            end_lines=end_lines)
    else:
        h_lines = ['DTASelect v2.1.12\n',
                   tabulated_args,
                   '\n',
                   'Locus	Sequence Count	Spectrum Count	Sequence Coverage	Length	MolWt	pI	Validation Status	'
                   'NSAF	EMPAI	Descriptive Name	HRedundancy	LRedundancy	MRedundancy\n',
                   'Unique	FileName	XCorr	DeltCN	Conf%	M+H+	CalcM+H+	PPM	TotalIntensity	SpR	Prob Score	'
                   'pI	IonProportion	Redundancy	Sequence	RetTime	PTMIndex	PTMIndex Protein List\n'
                   ]
        dta_filter_content = dtaselectfilter.to_dta_select_filter(
            version=dtaselectfilter.DtaSelectFilterVersion.V2_1_12_paser,
            h_lines=h_lines,
            dta_filter_results=filter_results,
            end_lines=end_lines)

    return alignment_figs, pin_df, StringIO(dta_filter_content)


if __name__ == '__main__':
    run()
