from __future__ import annotations

import math
from io import StringIO, TextIOWrapper
from serenipy import sqt, dtaselectfilter
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from .util import calculate_protein_coverage, get_unmodified_peptide, map_protein_to_peptides, map_peptide_to_specid


def convert_to_csv(sqt_content: TextIOWrapper | StringIO, filename=None) -> pd.DataFrame:
    """
    Convert SQT filee into a pandas df
    :param sqt_content: contents of the sqt file
    :param filename: name of the sqt file
    :return: the pandas Dataframe
    """

    _, _, s_lines = sqt.from_sqt(sqt_content)

    data = {'low_scan': [], 'high_scan': [], 'charge': [], 'process_time': [], 'server': [],
            'experimental_mass': [], 'total_ion_intensity': [], 'lowest_sp': [], 'number_matches': [],
            'experimental_ook0': [], 'experimental_mz': [], 'corrected_ook0': [],

            'xcorr_rank': [], 'sp_rank': [], 'calculated_mass': [], 'delta_cn': [], 'xcorr': [], 'sp': [],
            'matched_ions': [],
            'expected_ions': [], 'sequence': [], 'validation_status': [], 'predicted_ook0': [], 'tims_score': [],
            'tims_b_score_m2': [],
            'tims_b_score_best_m': [],

            'locuses': [], 'target': [], 'm_line': []}

    for s_line in s_lines:
        for i, m_line in enumerate(s_line.m_lines[:-1]):
            data['low_scan'].append(s_line.low_scan)
            data['high_scan'].append(s_line.high_scan)
            data['charge'].append(s_line.charge)
            data['process_time'].append(s_line.process_time)
            data['server'].append(s_line.server)
            data['experimental_mass'].append(s_line.experimental_mass)
            data['total_ion_intensity'].append(s_line.total_ion_intensity)
            data['lowest_sp'].append(s_line.lowest_sp)
            data['number_matches'].append(s_line.number_matches)
            data['experimental_ook0'].append(s_line.experimental_ook0)
            data['experimental_mz'].append(s_line.experimental_mz)
            data['corrected_ook0'].append(s_line.corrected_ook0)

            data['xcorr_rank'].append(m_line.xcorr_rank)
            data['sp_rank'].append(m_line.sp_rank)
            data['calculated_mass'].append(m_line.calculated_mass)
            data['delta_cn'].append(s_line.m_lines[i + 1].delta_cn)
            data['xcorr'].append(m_line.xcorr)
            data['sp'].append(m_line.sp)
            data['matched_ions'].append(m_line.matched_ions)
            data['expected_ions'].append(m_line.expected_ions)
            data['sequence'].append(m_line.sequence)
            data['validation_status'].append(m_line.validation_status)
            data['predicted_ook0'].append(m_line.predicted_ook0)
            data['tims_score'].append(m_line.tims_score)
            data['tims_b_score_m2'].append(m_line.tims_b_score_m2)
            data['tims_b_score_best_m'].append(m_line.tims_b_score_best_m)

            data['locuses'].append(' '.join([l_line.locus_name for l_line in m_line.l_lines]))

            data['target'].append(not all(['Reverse' in l_line.locus_name for l_line in m_line.l_lines]))
            data['m_line'].append(i)

    df = pd.DataFrame(data)
    df.index.name = 'id'
    if filename:
        df['file'] = filename
    return df


def convert_to_moka(sqt_df: pd.DataFrame):
    """
    Converts a sqt df to a mokapot input df (pin format)
    :param sqt_df: sqt df
    :return: mokapot input df (pin)
    """
    perc_df = pd.DataFrame()
    perc_df.index = sqt_df.index
    perc_df.index.name = 'SpecId'

    perc_df['Label'] = sqt_df['target']
    perc_df.loc[perc_df['Label'] == 0, 'Label'] = -1
    perc_df.loc[perc_df['Label'] == 1, 'Label'] = 1

    perc_df['ScanNr'] = sqt_df['low_scan']
    perc_df['ExpMass'] = sqt_df['experimental_mass']
    perc_df['CalcMass'] = sqt_df['calculated_mass']

    # perc_df['group'] = df['file']

    perc_df['abs_ppm'] = \
        abs((sqt_df['experimental_mass'] - sqt_df['calculated_mass']).divide(sqt_df['calculated_mass']) * 1_000_000)
    perc_df['abs_mass_diff'] = abs(sqt_df['experimental_mass'] - sqt_df['calculated_mass'])

    charges = pd.get_dummies(sqt_df['charge'], prefix='charge')
    perc_df = pd.concat([perc_df, charges], axis=1)

    perc_df['delta_cn'] = sqt_df['delta_cn']
    perc_df['xcorr'] = sqt_df['xcorr']
    perc_df['matched_ions'] = sqt_df['matched_ions']
    perc_df['expected_ions'] = sqt_df['expected_ions']
    perc_df['matched_ion_fraction'] = sqt_df['matched_ions'].divide(sqt_df['expected_ions'])

    perc_df['sp'] = sqt_df['sp']
    perc_df['sequence_length'] = [len(get_unmodified_peptide(peptide)) for peptide in sqt_df['sequence']]

    perc_df['xcorr_rank'] = sqt_df['xcorr_rank']
    perc_df['sp_rank'] = sqt_df['sp_rank']
    perc_df['tryp'] = [int(peptide[0] in ['KR']) + int(peptide[-3] in ['KR']) for peptide in sqt_df['sequence']]

    m_lines = pd.get_dummies(sqt_df['m_line'], prefix='m_line')
    perc_df = pd.concat([perc_df, m_lines], axis=1)

    if any(sqt_df['tims_score']):
        perc_df['tims_score'] = sqt_df['tims_score']
        perc_df['tims_score'] = perc_df['tims_score'].fillna(value=0)

    perc_df['Peptide'] = sqt_df['sequence']
    perc_df['Proteins'] = sqt_df['locuses']
    perc_df.reset_index(inplace=True)

    return perc_df


def get_filter_results_moka(sqt_df: pd.DataFrame, psm_results: pd.DataFrame, peptide_results: pd.DataFrame,
                            protein_results: pd.DataFrame, fasta_dict) -> list[dtaselectfilter.DTAFilterResult]:
    """
    Converts mokapot results into DTAFilterResult results
    :param sqt_df: sqt df
    :param psm_results: mokapot psm results
    :param peptide_results: mokapot peptide results
    :param protein_results: mokapot protein results
    :param fasta_dict: locus to protein sequence map
    :return: list of DTAFilterResult results
    """
    peptide_to_specid = map_peptide_to_specid(psm_results)
    protein_to_peptides = map_protein_to_peptides(protein_results, peptide_results)

    filter_results = []
    for protein, peptides in protein_to_peptides.items():
        protein_lines = []
        peptide_lines = []
        peptides = [peptide for peptide in peptides if peptide in peptide_to_specid]
        psm_ids = [psm_id for peptide in peptides for psm_id in peptide_to_specid[peptide]]
        sqt_psm_df = sqt_df.loc[psm_ids]
        for locus in [protein]:
            protein_sequence = fasta_dict[locus]['sequence']
            description_name = fasta_dict[locus]['description']
            sequence_count = len({(sequence, charge) for sequence, charge in sqt_psm_df[['sequence', 'charge']].values})
            spectrum_count = len(psm_ids)
            sequence_coverage = sum(calculate_protein_coverage(protein_sequence, peptides)) / len(
                protein_sequence) * 100
            length = len(protein_sequence)

            try:
                protein = ProteinAnalysis(protein_sequence)
                molWt = protein.molecular_weight()
                pi = protein.isoelectric_point()
            except ValueError:
                molWt = None
                pi = None

            protein_line = dtaselectfilter.ProteinLine(locus_name=locus,
                                                       sequence_count=sequence_count,
                                                       spectrum_count=spectrum_count,
                                                       sequence_coverage=sequence_coverage,
                                                       length=length,
                                                       molWt=molWt,
                                                       pi=pi,
                                                       validation_status='U',
                                                       nsaf=None,
                                                       empai=(math.pow(10, (sequence_coverage / 100))) - 1,
                                                       description_name=description_name,
                                                       h_redundancy=None,
                                                       l_redundancy=None,
                                                       m_redundancy=None)
            protein_lines.append(protein_line)

        for _, row in sqt_psm_df.iterrows():
            peptide = ProteinAnalysis(row['sequence'])

            file_path = row['file']
            low_scan = row['low_scan']
            high_scan = row['high_scan']
            charge = row['charge']

            file_name = f'{file_path}.{low_scan}.{high_scan}.{charge}'

            try:
                ppm = (row['experimental_mass'] - row['calculated_mass']) / row['calculated_mass'] * 1_000_000
            except ZeroDivisionError:
                ppm = None

            try:
                ion_proportion = row['matched_ions'] / row['expected_ions'] * 100
            except ZeroDivisionError:
                ion_proportion = None

            peptide_line = dtaselectfilter.PeptideLine(
                unique=None,
                file_name=file_name,
                x_corr=row['xcorr'],
                delta_cn=row['delta_cn'],
                conf=None,
                mass_plus_hydrogen=row['experimental_mass'],
                calc_mass_plus_hydrogen=row['calculated_mass'],
                ppm=ppm,
                total_intensity=row['total_ion_intensity'],
                spr=row['sp_rank'],
                ion_proportion=ion_proportion,
                redundancy=None,  # add in later
                sequence=row['sequence'],
                prob_score=row['sp'],
                pi=peptide.isoelectric_point(),
                measured_im_value=row['experimental_ook0'],
                predicted_im_value=row['predicted_ook0'],
                im_score=row['tims_score'],
                ret_time=None,
                ptm_index=None,
                ptm_index_protein_list=None,
                experimental_mz=row['experimental_mz'],
                corrected_1k0=row['corrected_ook0'],
                ion_mobility=None)

            peptide_lines.append(peptide_line)
        peptide_lines.sort(key=lambda x: x.low_scan)
        filter_results.append(dtaselectfilter.DTAFilterResult(protein_lines=protein_lines, peptide_lines=peptide_lines))
    return filter_results
