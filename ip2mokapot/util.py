import os
import re

import numpy as np


def map_locus_to_sequence_from_fasta(fasta_lines):
    locus_to_sequence_map = {}
    locus = None
    for line in fasta_lines:
        if line == "":
            continue
        elif line[0] == ">":  # new protein
            locus = line.rstrip().split(" ")[0].replace(">", "")
            description = " ".join(line.rstrip().split(" ")[1:])
            locus_to_sequence_map[locus] = {'sequence': "", 'description': description}
        else:  # protein sequence
            locus_to_sequence_map[locus]['sequence'] += line.rstrip()
    return locus_to_sequence_map


def get_unmodified_peptide(peptide_sequence: str) -> str:
    if peptide_sequence[2] == '.' and peptide_sequence[-2] == '.':
        peptide_sequence = peptide_sequence[2:-2]

    pattern = re.compile(r'[^A-Z]')
    return pattern.sub('', peptide_sequence)

def calculate_protein_coverage(protein, peptides):
    cov_arr = np.zeros(len(protein))
    for peptide in peptides:

        if peptide[2] == '.' and peptide[-2] == '.':
            peptide = peptide[2:-2]

        peptide = get_unmodified_peptide(peptide)

        for i in [i.start() for i in re.finditer(peptide, protein)]:
            cov_arr[i:i + len(peptide)] = 1
    return cov_arr


def xml_to_dict(xml_file):
    import xml.etree.ElementTree as ET

    tree = ET.parse(xml_file)
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


def run_perc(input_file, fasta_file, perc_path):
    decoy_flag = "Reverse_"

    target_psms = "target_psms.tsv"
    target_peptides = "target_peptides.tsv"
    target_proteins = "target_proteins.tsv"
    decoy_psms = "decoy_psms.tsv"
    decoy_peptides = "decoy_peptides.tsv"
    decoy_proteins = "decoy_proteins.tsv"

    fasta_arg_str = f"-f {fasta_file}"
    perc_command = f"{perc_path} -j {input_file} -P {decoy_flag} " \
                   f"{fasta_arg_str if fasta_file else '-A'} -A " \
                   f"-m {target_psms} -r {target_peptides} -l {target_proteins} " \
                   f"-M {decoy_psms} -B {decoy_peptides} -L {decoy_proteins}"

    print(perc_command)
    os.system(perc_command)

    return target_psms, target_peptides, target_proteins, decoy_psms, decoy_peptides, decoy_proteins


def map_protein_to_peptides(protein_results, peptide_results):
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
    peptide_to_specid = {}
    for specid, peptide in psm_results[['SpecId', 'Peptide']].values:
        peptide_to_specid.setdefault(peptide, set()).add(specid)
    return peptide_to_specid