import re
from collections import defaultdict
from xml.etree import ElementTree

import numpy as np
from mokapot.parsers.fasta import _group_proteins, _cleave
from mokapot.proteins import Proteins, LOGGER


def get_unmodified_peptide(peptide_sequence: str) -> str:
    """
    Converts a peptide sequence to an unmodified version. cleans n-term and c-term amino acids if they are present.
    :param peptide_sequence: peptides sequence
    :return: unmodified peptide sequence
    """
    if peptide_sequence[2] == '.' and peptide_sequence[-2] == '.':
        peptide_sequence = peptide_sequence[2:-2]

    pattern = re.compile(r'[^A-Z]')
    return pattern.sub('', peptide_sequence)


def calculate_protein_coverage(protein: str, peptides: list[str]) -> np.ndarray:
    """
    Calculates the protein sequence coverage
    :param protein: protein sequence
    :param peptides: list of peptides
    :return: np array which specifies where the protein is covered
    """
    cov_arr = np.zeros(len(protein))
    for peptide in peptides:

        if peptide[2] == '.' and peptide[-2] == '.':
            peptide = peptide[2:-2]

        peptide = get_unmodified_peptide(peptide)

        for i in [i.start() for i in re.finditer(peptide, protein)]:
            cov_arr[i:i + len(peptide)] = 1
    return cov_arr


def xml_to_dict(xml_file):
    """
    Converts a xml file to python dict
    :param xml_file: xml file
    :return: dict containing xml file values
    """

    tree = ElementTree.parse(xml_file)
    root = tree.getroot()
    data = {}

    def recurse_through_xml(root, data):
        for child in root:
            if len(child) == 0:
                if child.tag in data:
                    if not isinstance(data[child.tag], list):
                        data[child.tag] = [data[child.tag], child.text]
                    else:
                        data[child.tag].append(child.text)
                else:
                    data[child.tag] = child.text
            else:
                if child.tag in data:
                    if not isinstance(data[child.tag], list):
                        data[child.tag] = [data[child.tag], {}]
                        recurse_through_xml(child, data[child.tag][-1])
                    else:
                        data[child.tag].append({})
                        recurse_through_xml(child, data[child.tag][-1])
                else:
                    data[child.tag] = {}
                    recurse_through_xml(child, data[child.tag])

    recurse_through_xml(root, data)
    return data


def map_protein_to_peptides(protein_results, peptide_results):
    """
    maps protein locus to peptide sequence
    :param protein_results: mokapot protein results df
    :param peptide_results: mokapot peptide results df
    :return: map of protein locus to peptide sequence
    """
    protein_to_peptide_map = {}
    for proteins in protein_results['mokapot protein group']:
        proteins = proteins.split(', ')
        for protein in proteins:
            protein_to_peptide_map.setdefault(protein, set())

    for peptide, proteins in peptide_results[['Peptide', 'Proteins']].values:
        for protein in proteins.split(' '):
            if protein in protein_to_peptide_map:
                protein_to_peptide_map[protein].add(peptide)

    return protein_to_peptide_map


def map_peptide_to_specid(psm_results):
    """
    maps peptide sequence to specid
    :param psm_results: mokapot psm results df
    :return: peptide to specid map
    """
    peptide_to_specid = {}
    for specid, peptide in psm_results[['SpecId', 'Peptide']].values:
        peptide_to_specid.setdefault(peptide, set()).add(specid)
    return peptide_to_specid


def strip_modifications(peptide_sequence: str) -> str:
    """
    Removes any non-amino-acid characters from the given peptide sequence.

    Args:
        peptide_sequence: The peptide sequence to be stripped of modifications.

    Returns:
        The peptide sequence with all non-amino-acid characters removed.
    """
    return ''.join([c for c in peptide_sequence if c in 'ACDEFGHIKLMNPQRSTVWY'])
