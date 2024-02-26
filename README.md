# Syn-CpG-Spacer

Syn-CpG-Spacer is a Progressive Web App (PWA) for biomedical scientists written in Python using Panel, Bokeh, Biopython libraries. It allows for synonymous recoding of genetic sequences to increase the frequency of CpG dinucleotides by setting constraints on their spacing. The primary usecase are experiments with attenuation of viruses.

The software changes codons along a sequence to synonymous alternatives that form CpG dinucleotides according to the user's settings. This can be done at codon positions 1-2, 2-3 and 3-1 (split over two subsequent codons).

Using Panel's Pyodide integration, the app is hosted on GitHub Pages in this repository and is available on the following address:

https://oleksulkowski.github.io/Syn-CpG-Spacer/app/

## Installation

In browsers such as Chrome, Safari or Edge, it is possible to install the app onto your machine for offline use by clicking the browser prompt after opening the link. An installed app will download and apply updates automatically when they become available.

## Usage

The software allows the user to load their own FASTA sequence or to use a pre-loaded sample sequence (part of <a href="https://www.ncbi.nlm.nih.gov/nucleotide/MN685337.1">HIV-1 Gag</a>). The user can then either set a minimum gap between newly added CpG's or set a desired <ins>average</ins> gap between CpG's. With the latter option, the software will find a minimum gap that will result in as close a possible average gap to the user's setting using a binary search algorithm.

The program allows protecting a set number of initial and final nucleotides from changes, which might be biologically relevant. As increasing the CpG content can decrease the frequency of A in a sequence, the user can also decide to make the remaining sequence synonymously A-rich after CpG's have been added.

Every new recoded sequence requires input of a unique ID. The sequences are displayed on an interactive alignment view that highlights CpG dinucleotides. A table shows statistical data. The user can adjust the settings and compare the sequences. When finished, the user can download the outputs as a FASTA file.

## Algorithm outline

1. The user configures the minimum CpG gap, protected terminal nucleotide length and chooses whether to make the sequence A-rich after adding CpG's.
    - If the user sets a target average gap, a binary search algorithm will perform the steps below to find a minimum CpG gap that results in the closest average CpG gap to the desired one.
2. Codon instances are generated for every codon along the sequence. It is checked whether the codon already contains a CpG or forms a split CpG with the next codon.
3. It is determined which codons can potentially be transformed into CpG-forming alternatives based on their position in the sequence. The criterium is being at least the minimum CpG gap away from existing CpG's.
4. The initial and final number of nucleotides are protected against changes, if specified by the user.
5. Codons are mutated to synonymous CpG-forming alternatives along the sequence. Minimum CpG gap between newly added CpG's is ensured.
6. The sequences' synonymity is checked, along with the preservation of terminal signals and adherence to the minimum gap settings.
7. If the A-enrichment option is selected, the rest of the sequence is synonymously recoded into more A-rich codons, without impacting CpG's.
8. The same checks as those described in step 6 are performed.


## Development

Use the `environment.yml` file to create an environment with all the dependencies:

```
conda env create -f environment.yml
conda activate Syn-CpG-Spacer
```

As per Panel <a href="https://panel.holoviz.org/how_to/wasm/">documentation</a>, develop locally in `index.py` using
```
panel serve index.py --autoreload
```

After making changes, convert `index.py` to the Pyodide PWA:
```
panel convert index.py --to pyodide-worker --out docs/app --title Syn-CpG-Spacer --pwa
```


You can run the Pyodide app locally on http://localhost:8000/docs/app by using
```
python3 -m http.server
```

## Tests

Syn-CpG-Spacer uses Pytest for checking if code changes introduced errors into the recoding algorithm by comparing the new output to a set of validated sequences. This is hooked up to Github Actions CI. Run the tests using

```
pytest
```

Within the app, each algorithm run is checked to ensure correct application of user-defined variables.


## Community contributions

Please use the <a href="https://github.com/oleksulkowski/synrecoder/issues">issues tab</a> for bug reports and feature requests.

## Acknowledgements

The Bokeh sequence viewer is based on <a href="https://dmnfarrell.github.io/bioinformatics/bokeh-sequence-aligner">code</a> by Damien Farrell (<a href="https://github.com/dmnfarrell">@dmnfarrell</a>).