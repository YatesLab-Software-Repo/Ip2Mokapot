import re
from collections import defaultdict
from xml.etree import ElementTree
import multiprocessing as mp

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


def _parse_fasta_files(fasta_files):
    """Read a fasta file and divide into proteins

    Parameters
    ----------
    fasta_files : str or list of str
        One or more FASTA files.

    Returns
    -------
    proteins : list of str
        The raw protein headers and sequences.
    """
    fasta_files = fasta_files
    fasta = []
    for fa in fasta_files:
        fasta.append(fa.read())

    return "\n".join(fasta)[1:].split("\n>")


def _parse_protein(raw_protein):
    """Parse the raw string for a protein.

    Parameters
    ----------
    raw_protein : str
        The raw protein string.

    Returns
    -------
    header : str
        The protein name.
    sequence : str
        The protein sequence.
    """
    entry = raw_protein.splitlines()
    prot = entry[0].split(" ")[0]
    desc = "".join(entry[0].split(" ")[1:])
    if len(entry) == 1:
        return prot, "", desc

    seq = "".join(entry[1:])
    return prot, seq, desc


def _cleavage_sites(sequence, enzyme_regex, enzyme_term):
    """Find the cleavage sites in a sequence.

    Parameters
    ----------
    sequence : str
        A protein sequence to digest.
    enzyme_regex : str or compiled regex
        A regular expression defining the enzyme specificity.

    Returns
    -------
    sites : list of int
        The cleavage sites in the sequence.
    """
    if isinstance(enzyme_regex, str):
        enzyme_regex = re.compile(enzyme_regex)

    # Find the cleavage sites
    sites = (
        [0]
        + [m.end() if enzyme_term else m.start() for m in enzyme_regex.finditer(sequence)]
        + [len(sequence)]
    )
    return sites

def digest(
    sequence,
    enzyme_regex="[KR]",
    missed_cleavages=0,
    clip_nterm_methionine=False,
    min_length=6,
    max_length=50,
    semi=False,
    enzyme_term=True
):
    """
    Digest a protein sequence into its constituent peptides.

    Parameters
    ----------
    sequence : str
        A protein sequence to digest.
    enzyme_regex : str or compiled regex, optional
        A regular expression defining the enzyme specificity. The end of the
        match should indicate the cleavage site.
    missed_cleavages : int, optional
        The maximum number of allowed missed cleavages.
    clip_nterm_methionine : bool, optional
        Remove methionine residues that occur at the protein N-terminus.
    min_length : int, optional
        The minimum peptide length.
    max_length : int, optional
        The maximum peptide length.
    semi : bool
        Allow semi-enzymatic cleavage.

    Returns
    -------
    peptides : set of str
        The peptides resulting from the digested sequence.
    """
    sites = _cleavage_sites(sequence, enzyme_regex, enzyme_term)

    peptides = _cleave(
        sequence=sequence,
        sites=sites,
        missed_cleavages=missed_cleavages,
        min_length=min_length,
        max_length=max_length,
        semi=semi,
        clip_nterm_met=clip_nterm_methionine,
    )

    return peptides


def process_protein(seq, prot, enzyme_regex, missed_cleavages, min_length, max_length, semi, clip_nterm_methionine, enzyme_term):
    peps = digest(
        seq,
        enzyme_regex=enzyme_regex,
        missed_cleavages=missed_cleavages,
        min_length=min_length,
        max_length=max_length,
        semi=semi,
        clip_nterm_methionine=clip_nterm_methionine,
        enzyme_term=enzyme_term
    )

    return (peps, prot)

def read_fasta(
        fasta,
        enzyme="[KR]",
        missed_cleavages=2,
        clip_nterm_methionine=False,
        min_length=6,
        max_length=50,
        semi=False,
        decoy_prefix="decoy_",
        enzyme_term=True
):
    """Parse a FASTA file, storing a mapping of peptides and proteins.

    Protein sequence information from the FASTA file is required to compute
    protein-level confidence estimates using the picked-protein approach.
    Decoys proteins must be included and must be of the have a description in
    format of `<prefix><protein ID>` for valid confidence estimates to be
    calculated.

    If you need to generate an appropriate FASTA file with decoy sequences for
    your database search, see :py:func:`mokapot.make_decoys()`.

    Importantly, the parameters below should match the conditions in which the
    PSMs were assigned as closely as possible. Enzyme specificity is provided
    using a regular expression. A table of common enzymes can be found here in
    the mokapot `cookbook
    <file:///Users/wfondrie/packages/mokapot/docs/build/html/cookbook.html#enzyme-regular-expressions>`_.

    Parameters
    ----------
    fasta : list of TextIOWrapper's and StringIO's
        The FASTA file(s) used for assigning the PSMs
    decoy_prefix : str, optional
        The prefix used to indicate a decoy protein in the description
        lines of the FASTA file.
    enzyme : str or compiled regex, optional
        A regular expression defining the enzyme specificity was used when
        assigning PSMs. The cleavage site is interpreted as the end of the
        match. The default is trypsin, without proline suppression: "[KR]".
    missed_cleavages : int, optional
        The allowed number of missed cleavages.
    clip_nterm_methionine : bool, optional
        Remove methionine residues that occur at the protein N-terminus.
    min_length : int, optional
        The minimum peptide length to consider.
    max_length : int, optional
        The maximum peptide length to consider.
    semi : bool, optional
        Was a semi-enzymatic digest used to assign PSMs? If :code:`True`, the
        protein database will likely contain many shared peptides and yield
        unhelpful protein-level confidence estimates.

    Returns
    -------
    Proteins object
        The parsed proteins as a :py:class:`~mokapot.proteins.Proteins`
        object.

    """
    if isinstance(enzyme, str):
        enzyme_regex = re.compile(enzyme)
    else:
        enzyme_regex = enzyme

    # Read in the fasta files
    LOGGER.info("Parsing FASTA files and digesting proteins...")
    # Build the initial mapping
    proteins = {}
    peptides = defaultdict(set)

    pool = mp.Pool(processes=mp.cpu_count())
    results = [pool.apply_async(process_protein, args=(seq,prot, enzyme_regex, missed_cleavages, min_length, max_length, semi, clip_nterm_methionine, enzyme_term)) for prot, seq, desc in fasta]

    pool.close()
    pool.join()

    for result in results:
        peps, prot = result.get()
        if peps:
            proteins[prot] = peps
            for pep in peps:
                peptides[pep].add(prot)

    print(proteins)
    total_prots = len(fasta)
    LOGGER.info("  - Parsed and digested %i proteins.", total_prots)
    LOGGER.info("  - %i had no peptides.", len(fasta) - len(proteins))
    LOGGER.info("  - Retained %i proteins.", len(proteins))
    del fasta

    # Sort proteins by number of peptides:
    proteins = {
        k: v for k, v in sorted(proteins.items(), key=lambda i: len(i[1]))
    }

    LOGGER.info("Matching target to decoy proteins...")
    # Build the decoy map:
    decoy_map = {}
    no_decoys = 0
    has_decoys = False
    has_targets = False
    for prot_name in proteins:
        if not prot_name.startswith(decoy_prefix):
            has_targets = True
            decoy = decoy_prefix + prot_name
            decoy_map[prot_name] = decoy
            if decoy in proteins.keys():
                has_decoys = True
            else:
                no_decoys += 1

    if not has_targets:
        raise ValueError("Only decoy proteins were found in the FASTA file.")

    if no_decoys and no_decoys < len(decoy_map):
        LOGGER.warning(
            "Found %i target proteins without matching decoys.", no_decoys
        )

    LOGGER.info("Building protein groups...")
    # Group Proteins
    num_before_group = len(proteins)
    proteins, peptides = _group_proteins(proteins, peptides)
    LOGGER.info(
        "\t- Aggregated %i proteins into %i protein groups.",
        num_before_group,
        len(proteins),
    )

    if not has_decoys:
        LOGGER.info("No decoy sequences were found in the FASTA file.")
        LOGGER.info(
            "  - Creating decoy protein groups that mirror the target "
            "proteins."
        )

    # unique peptides:
    LOGGER.info("Discarding shared peptides...")
    shared_peptides = {}
    unique_peptides = {}
    for pep, prots in peptides.items():
        if len(prots) == 1:
            unique_peptides[pep] = next(iter(prots))
        else:
            shared_peptides[pep] = "; ".join(prots)

    total_proteins = len(set(p for p in unique_peptides.values()))

    LOGGER.info(
        "  - Discarded %i peptides and %i proteins groups.",
        len(peptides) - len(unique_peptides),
        len(proteins) - total_proteins,
    )

    LOGGER.info(
        "  - Retained %i peptides from %i protein groups.",
        len(unique_peptides),
        total_proteins,
    )

    parsed = Proteins(
        decoy_prefix=decoy_prefix,
        peptide_map=unique_peptides,
        shared_peptides=shared_peptides,
        protein_map=decoy_map,
        has_decoys=has_decoys,
    )

    return parsed

def strip_modifications(peptide_sequence: str) -> str:
    """
    Removes any non-amino-acid characters from the given peptide sequence.

    Args:
        peptide_sequence: The peptide sequence to be stripped of modifications.

    Returns:
        The peptide sequence with all non-amino-acid characters removed.
    """
    return ''.join([c for c in peptide_sequence if c in 'ACDEFGHIKLMNPQRSTVWY'])