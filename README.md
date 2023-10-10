# Syn-CpG-Spacer

Syn-CpG-Spacer is a Progressive Web App (PWA) for biomedical scientists written in Python using Panel, Bokeh, Biopython libraries. It allows for synonymous recoding of genetic sequences to increase the frequency of CpG dinucleotides by setting constraints on their spacing. The primary usecase are experiments with attenuation of viruses.

The software changes codons along a sequence to synonymous alternatives that form CpG dinucleotides according to user's settings. This can be done at codon positions 1-2, 2-3 and 3-1 (split over two subsequent codons).

Thanks to Panel's Pyodide integration, the app is hosted on GitHub Pages in this repository and is available at the following address:

https://oleksulkowski.github.io/app

## Installation

In browsers such as Chrome, Safari or Edge, it is possible to install the app onto your machine for offline use by clicking the browser prompt after opening the link. An installed app will download and apply updates automatically when they become available.

## Usage

The software allows the user to load their own FASTA sequence or to use a pre-loaded sample sequence (part of <a href="https://www.ncbi.nlm.nih.gov/nucleotide/MN685337.1">HIV-1 Gag</a>). The user can then either set a minimum gap between newly added CpG's or set a desired average gap between CpG's. With the latter option, the software will find a minimum gap that will result in as close a possible average gap to the user's setting using a binary search algorithm.

Due to the fact that in some viral genomes, mutations in terminal regions can interfere with packaging, the program allows protecting a set number of initial and final nucleotides from changes. As increasing the CpG content can decrease the frequency of A in a sequence, the user can also decide to make the remaining sequence synonymously more A-rich after CpG's have been added.

Every new recoded sequence requires input of a unique ID. The sequences are displayed on an interactive alignment view that highlights CpG dinucleotides. A table shows statistical data. The user can adjust the settings and compare the sequences. When finished, the user can download the outputs as a FASTA file.

## Algorithm outline

1. The user configures the minimum CpG gap, protected terminal nucleotide length and chooses whether to make the sequence A-rich after adding CpG's.
    - If the user sets a target average gap, a binary search algorithm will perform the steps below to find a minimum CpG gap that results in the closest average CpG gap to the desired one.
2. Codon instances are generated for every codon along the sequence. It is checked whether the codon already contains a CpG or forms a split CpG with the next codon.
3. It is determined which codons can potentially be transformed into CpG-forming alternatives based on their position in the sequence. The criterium is being at least the minimum CpG gap away from existing CpG's.
4. The initial and final number of nucleotides are protected against changes, if specified by the user.
5. Codons are mutated to synonymous CpG-forming alternatives along the sequence. Minimum CpG gap between newly added CpG's is ensured.
6. The sequences' synonymity is checked, along with the preservation of packaging signals and adherence to the minimum gap settings.
7. If the A-enrichment option is selected, the rest of the sequence is synonymously recoded into more A-rich codons, without impacting CpG's.
8. The same checks as those described in step 6 are performed.

## Community guidelines


## Tests

