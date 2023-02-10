# Ip2Mokapot (MokaFilter)

Ip2Mokapot is a software program that serves as a container for the Mokapot analysis pipeline. It allows users to seamlessly integrate their IP2 input files into the Mokapot analysis process, and output results in an IP2 compatible format. This program provides a convenient and easy-to-use interface for users who would like to take advantage of the powerful analysis capabilities of Mokapot without having to modify their existing workflows. It's important to note that Ip2Mokapot does include the Mokapot package and is simply a tool that facilitates the use of Mokapot in an IP2 environment.

## Citing Mokapot

> Fondrie W. E. & Noble W. S. mokapot: Fast and Flexible Semisupervised
> Learning for Peptide Detection. J Proteome Res (2021) doi:
> 10.1021/acs.jproteome.0c01010. PMID: 33596079.
> [Link](https://doi.org/10.1021/acs.jproteome.0c01010)


## How to install:

>$pip install git+https://github.com/YatesLab-Software-Repo/Ip2MokaPot

## How to run (IP2):

>$mokafilter --sqts data\test.sqt --fastas data\test.fasta --search_xml data\search.xml --out data\DTASelect-filter.txt --dta_params data\DTASelect.params

## How to run (PaSER):

>$mokafilter --sqts data\test.sqt --fastas data\test.fasta --search_xml data\search.xml --out data\DTASelect-filter.txt --dta_params data\DTASelect.params --timscore True

## See arguments:

>$mokafilter -h

## Run local Streamlit app

>$streamlit run home.py
