import pandas as pd
import numpy as np
from io import StringIO, BytesIO
import math
import copy


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d, SingleIntervalTicker
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot


import panel as pn
import panel.widgets as pnw

pn.extension("tabulator", "notifications", design="fast", notifications=True)
template = pn.template.FastListTemplate(
    title="Syn-CpG-Spacer",
    theme_toggle=False,
    accent="#A01346",
    collapsed_sidebar=True,
    main_max_width="80vw",
)

notifications = pn.state.notifications

# Taken from github.com/Mmark94/protein_recoding. For STOP codons, X is used as symbol.
dna_to_pro = {
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
    "TGC": "C",
    "TGT": "C",
    "GAC": "D",
    "GAT": "D",
    "GAA": "E",
    "GAG": "E",
    "TTC": "F",
    "TTT": "F",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGT": "G",
    "CAC": "H",
    "CAT": "H",
    "ATA": "I",
    "ATC": "I",
    "ATT": "I",
    "AAA": "K",
    "AAG": "K",
    "CTA": "L",
    "CTC": "L",
    "CTG": "L",
    "CTT": "L",
    "TTA": "L",
    "TTG": "L",
    "ATG": "M",
    "AAC": "N",
    "AAT": "N",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCT": "P",
    "CAA": "Q",
    "CAG": "Q",
    "AGA": "R",
    "AGG": "R",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGT": "R",
    "AGC": "S",
    "AGT": "S",
    "TCA": "S",
    "TCC": "S",
    "TCG": "S",
    "TCT": "S",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACT": "T",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "TGG": "W",
    "TAA": "X",
    "TAG": "X",
    "TGA": "X",
    "TAC": "Y",
    "TAT": "Y",
}

DIC = {
    "ACA": "ACG",
    "ACC": "ACG",
    "ACT": "ACG",
    "CCA": "CCG",
    "CCC": "CCG",
    "CCT": "CCG",
    "AGA": "CGA",
    "AGG": "CGA",
    "GCA": "GCG",
    "GCC": "GCG",
    "GCT": "GCG",
    "AGC": "TCG",
    "AGT": "TCG",
    "TCA": "TCG",
    "TCC": "TCG",
    "TCT": "TCG",
}

DIC_for_split = {
    "AAT": "AAC",
    "ACA": "ACC",
    "ACG": "ACC",
    "ACT": "ACC",
    "AGT": "AGC",
    "ATA": "ATC",
    "ATT": "ATC",
    "CAT": "CAC",
    "CCA": "CCC",
    "CCT": "CCC",
    "AGA": "CGC",
    "AGG": "CGC",
    "CGA": "CGC",
    "CGG": "CGC",
    "CGT": "CGC",
    "CTA": "CTC",
    "CTG": "CTC",
    "CTT": "CTC",
    "TTA": "CTC",
    "TTG": "CTC",
    "GAT": "GAC",
    "GGA": "GGC",
    "GGG": "GGC",
    "GGT": "GGC",
    "GTA": "GTC",
    "GTG": "GTC",
    "GTT": "GTC",
    "TAT": "TAC",
    "TCA": "TCC",
    "TCT": "TCC",
    "TGT": "TGC",
    "TTT": "TTC",
}


DIC_for_A_rich = {
    "ACC": "ACA",
    "ACT": "ACA",
    "AGG": "AGA",
    "ATC": "ATA",
    "ATT": "ATA",
    "CAG": "CAA",
    "CCC": "CCA",
    "CCT": "CCA",
    "CGC": "CGA",
    "CGG": "CGA",
    "CGT": "CGA",
    "CTC": "CTA",
    "CTG": "CTA",
    "CTT": "CTA",
    "GAG": "GAA",
    "GCC": "GCA",
    "GCT": "GCA",
    "GGG": "GGA",
    "GGC": "GGA",
    "GGT": "GGA",
    "GTC": "GTA",
    "GTG": "GTA",
    "GTT": "GTA",
    "TCC": "TCA",
    "TCT": "TCA",
    "TTG": "TTA",
    # STOP codons are not mutated
}


class Codon:
    """
    Represents a codon, which is a sequence of three nucleotides that together form a unit of genetic code in a DNA or RNA molecule.

    Attributes:
        position (int): The position of the codon within the Gene object, indicating the index of the first nucleotide of the codon.
        sequence (str): The sequence of nucleotides in the codon.
        next_codon (Codon, optional): Reference to the next codon object in the sequence, if applicable.
        has_cg (bool): Indicates whether the codon contains a CpG dinucleotide.
        protected (bool): Flag indicating whether the codon should not be mutated.
        cg_split (bool): Indicates if a CpG dinucleotide is split over two subsequent codons.
        mutable (bool, optional): Indicates if the codon can be mutated to contain a CpG dinucleotide.
        potential_new_seq (str, optional): Holds a potential new sequence for the codon that might be considered for mutation, to form a CpG dinucleotide.
        mutated (bool): Tracks whether the codon has been mutated.

    Methods:
        __init__(self, position, sequence, next_codon=None, has_cg=False, protected=False, cg_split=False, mutable=None, potential_new_seq=None, mutated=False):
            Initializes a Codon object with the provided attributes.

        get_cg_position(self, sequence):
            Determines the position of a CpG dinucleotide within the codon or across adjacent codons, relative to the 'C'.

        get_mutation_positions(self):
            Identifies positions within the codon where mutations are possible based on the potential new sequence.

        translate_codon(self):
            Translates the codon sequence into its corresponding amino acid symbol according to the genetic code.

        check_codon_viability(self):
            Checks if the codon sequence is viable for translation based on its sequence.
    """

    def __init__(
        self,
        position,
        sequence,
        next_codon=None,
        has_cg=False,
        protected=False,
        cg_split=False,
        mutable=None,
        potential_new_seq=None,
        mutated=False,
    ):
        """
        Initializes a Codon object with the provided attributes.

        Parameters:
            position (int): The position of the codon within the Gene object.
            sequence (str): The nucleotide sequence of the codon.
            next_codon (Codon, optional): The next codon in the sequence.
            has_cg (bool): True if the codon contains a CpG dinucleotide, else False.
            protected (bool): True if the codon is should not be mutated, else False.
            cg_split (bool): True if a CpG dinucleotide is split between this codon and the next, else False.
            mutable (bool, optional): True if the codon can potentially be mutated to contain a CpG, else False.
            potential_new_seq (str, optional): A potential new sequence for mutation, forming a CpG dinucleotide.
            mutated (bool): True if the codon has been mutated, else False.
        """

        self.position = position  # Position of the codon, not the CG!
        self.sequence = sequence
        self.next_codon = next_codon
        self.has_cg = has_cg
        self.protected = protected
        self.cg_split = cg_split
        self.mutable = mutable
        self.potential_new_seq = potential_new_seq
        self.mutated = mutated

        self.check_codon_viability()

    def get_cg_position(self, sequence):
        """
        Determines the position of a CpG dinucleotide within the codon or across adjacent codons, relative to the 'C'.

        Parameters:
            sequence (str): The nucleotide sequence of the codon to check for CpG presence.

        Returns:
            int or None: The position of the 'C' in the CpG dinucleotide, if present; otherwise, None.
        """

        if sequence is not None:
            # Check for a CG pair within this codon
            for i in range(3):
                if sequence[i] == "C" and i < 2 and sequence[i + 1] == "G":
                    return self.position + i

            # Check for a CG pair split across two codons
            if (
                self.next_codon is not None
                and sequence[2] == "C"
                and self.next_codon[0] == "G"
            ):
                return self.position + 2

        # No CG pair was found
        return

    def get_mutation_positions(self):
        """
        Identifies positions within the codon where mutations are possible based on the potential new sequence.

        Returns:
            list of int: A list of positions within the codon where mutations differ from the original sequence.
        """

        diff_positions = []

        for i in range(3):
            if self.sequence[i] != self.potential_new_seq[i]:
                diff_positions.append(self.position + i)
        return diff_positions

    def translate_codon(self):
        """
        Translates the codon sequence into its corresponding amino acid symbol according to the genetic code.

        Returns:
            str or None: The single-letter symbol of the amino acid if the codon is valid; otherwise, None.
        """

        if self.sequence in dna_to_pro:
            aa_symbol = dna_to_pro[self.sequence]
            return aa_symbol

    # Checks if the codon is viable. If not, this function will raise an error
    def check_codon_viability(self):
        """
        Checks if the codon sequence is viable for translation based on its sequence.

        Raises:
            Exception: If the codon sequence cannot be translated into a known amino acid.
        """

        aa_symbol = self.translate_codon()
        if aa_symbol is None:
            raise Exception(
                "Your sequence does not start in-frame \
                or some codons are not in the codon table."
            )


class Gene:
    """
    Represents a sequence of nucleotides in DNA or RNA, as loaded by the user.

    Attributes:
        original_sequence (str): The nucleotide sequence of the gene as originally provided.
            This serves as the baseline genetic material upon which mutations and analyses are performed.

        new_sequence (str): Holds the modified nucleotide sequence following mutation operations.
            This attribute is updated dynamically as mutations are applied to the original_sequence.
            It allows for comparison between the original and mutated gene sequences to assess the effects of mutations.

        sequence_length (int): The length of the original_sequence. This attribute is useful for
            operations that require knowledge of the gene's length, such as iterating over the sequence,
            analyzing codon distribution, and applying mutations within specified regions.

        packaging_signal_length_beginning (int): Specifies the length of a protected region at the
            beginning of the sequence that should not be mutated. This is important for preserving
            essential elements such as promoters or regulatory sequences that are critical for gene expression.

        packaging_signal_length_end (int): Specifies the length of a protected region at the end of
            the sequence that should not be mutated. Like the beginning signal, this ensures that
            essential regulatory or terminating sequences are maintained intact.

        gap_method (int): Indicates the method used to determine the spacing or gaps between
            CpG sites. This parameter guides the mutation strategies to either enforce a minimum gap length or
            adjust CpG spacing to approximate a target average gap. The gap_method can be specified as an integer:
            - 1: Indicates that a minimum gap length between CpG sites is specified.
            - 2: Indicates that a target average gap length between CpG sites is specified.

        minimum_CpG_gap (int): The minimum number of nucleotides that must separate added CpG sites in the
            mutated sequence.

        desired_CpG_gap (int, optional): An optional target average gap length between CpG sites.
            If specified, mutation strategies may aim to adjust CpG spacing to approximate this value.

        original_cg_positions (list): A list of positions within the original_sequence where CpG sites
            are located. This information is key to understanding the original CpG distribution.

        mutable_positions (list): A list of positions within the gene sequence deemed suitable for
            mutation, based on criteria such as CpG spacing, packaging signals, and other factors.
            This list is used to direct where mutations can be applied without compromising essential
            gene functions or regulatory elements.

        original_codons (list): A list of Codon objects representing the codons in the original_sequence.
            This attribute facilitates detailed analysis and manipulation at the codon level, allowing for
            targeted mutations and the study of codon usage patterns.

        current_codons (list): A list of Codon objects representing the codons in the sequence after mutations.
            It mirrors the original_codons and reflects changes made during mutation processes, serving as a
            basis for analyzing the mutated sequence's codon composition and distribution. Used to create new_sequence.

        original_average_gap (float): Represents the average gap length between CpG sites in the original_sequence.

    Methods:
        __init__(self, original_sequence, gap_method, packaging_signal_length_beginning=0, packaging_signal_length_end=0, minimum_CpG_gap=12, desired_CpG_gap=None):
            Initializes the Gene object with the provided sequence and configuration parameters.

        check_gap_method(self):
            Determines the gap method based on user configuration and initiates mutation or gap optimalization processes.

        analyze_codons(self):
            Analyzes the provided nucleotide sequence to identify codons and CpG sites, creating Codon objects for each.

        enforce_packaging_signal(self, codon):
            Ensures a given codon mutation respects the packaging signal constraints.

        load_new_sequence(self):
            Constructs the new nucleotide sequence after mutations have been applied.

        determine_changeable_CpG(self):
            Identifies codons that can be synonymously mutated to increase CpG sites while respecting gap constraints.

        translate_sequence(self, sequence):
            Translates a nucleotide sequence into its corresponding amino acid sequence.

        first_difference(self, str1, str2):
            Identifies the first difference between two amino acid sequences. Used for debugging.

        check_synonymity(self):
            Ensures that mutations applied to the nucleotide sequence do not alter the encoded amino acid sequence.

        check_minimum_gaps(self):
            Validates that the minimum gap constraint between CpG sites is respected in the mutated sequence.

        check_packaging_signals(self):
            Ensures that mutations do not alter the designated packaging signal regions of the sequence.

        mutate_CpG(self):
            Applies mutations to increase CpG sites while ensuring synonymity and respecting gap and packaging constraints.

        mutate_A_rich(self):
            Mutates the remaining sequence synonymously to increase the prevalence of adenine nucleotides, avoiding packaging signals and CpG sites.

        find_desired_gap(self, desired_gap):
            Find the closest possible gap setting to achieve a desired average gap between CpG sites.

        count_CpGs(self, sequence):
            Counts the number of CpG sites within the gene sequence.

        calculate_average_gap(self, mode="current", st_dev=False):
            Calculates the average length of gaps between CpG sites in the gene sequence.

        calculate_CpG_abundance_change(self):
            Calculates the percentage change in CpG abundance between the original gene sequence and the mutated gene sequence.

        calculate_A_abundance_change(self):
            Calculates the percentage change in the abundance of adenine between the original and new gene sequences.
    """

    def __init__(
        self,
        original_sequence,
        gap_method,
        packaging_signal_length_beginning=0,
        packaging_signal_length_end=0,
        minimum_CpG_gap=12,
        desired_CpG_gap=None,
    ):
        """
        Initialize a Gene object with specific configuration parameters.

        Parameters:
            original_sequence (str): The original nucleotide sequence of the gene.
            gap_method (int): The method to use for determining new CpG gaps for mutating the sequence.
            packaging_signal_length_beginning (int, optional): Length of the packaging signal at the beginning of the sequence. Defaults to 0.
            packaging_signal_length_end (int, optional): Length of the packaging signal at the end of the sequence. Defaults to 0.
            minimum_CpG_gap (int, optional): The minimum gap between CpG sites in the sequence. Defaults to 12.
            desired_CpG_gap (int, optional): The desired gap between CpG sites, used for adjusting the sequence accordingly.
        """

        self.original_sequence = original_sequence
        self.new_sequence = None

        self.sequence_length = len(original_sequence)
        self.packaging_signal_length_beginning = packaging_signal_length_beginning
        self.packaging_signal_length_end = packaging_signal_length_end
        self.gap_method = gap_method
        self.minimum_CpG_gap = minimum_CpG_gap
        self.desired_CpG_gap = desired_CpG_gap

        self.original_cg_positions = []
        self.mutable_positions = []

        self.original_codons = self.analyze_codons()
        self.current_codons = copy.deepcopy(self.original_codons)
        self.current_cg_positions = self.original_cg_positions.copy()
        self.original_average_gap = self.calculate_average_gap("original")

        self.check_gap_method()

    @property
    def minimum_CpG_gap_extra(self):
        return self.minimum_CpG_gap + 1

    def check_gap_method(self):
        """
        Determine and apply the gap method based on the user's configuration.

        This method decides between maintaining a minimum CpG gap or adjusting to a desired mean CpG gap, and initiates the corresponding process.
        """

        if (
            self.gap_method is None
        ):  # This is required for checking the sequence after it is loaded and before the user has mutated it
            return
        elif self.gap_method == 1 and self.desired_CpG_gap is None:
            self.determine_changeable_CpG()
            self.mutate_CpG()
        elif self.gap_method == 2 and self.minimum_CpG_gap is None:
            self.closest_gap = self.find_desired_gap(self.desired_CpG_gap)

    def analyze_codons(self):
        """
        Analyze the nucleotide sequence to identify codons and CpG sites.

        This method creates Codon objects for each codon in the sequence and identifies original CpG positions.
        """

        codons = []
        for i in range(0, self.sequence_length, 3):
            codon_sequence = self.original_sequence[i : i + 3]
            has_cg = "CG" in codon_sequence
            cg_split = None
            if i + 3 < len(self.original_sequence):
                next_codon = self.original_sequence[i + 3 : i + 6]
                if "CG" in (codon_sequence[-1] + next_codon[0]):
                    has_cg = True
                    cg_split = True
            else:
                has_cg = "CG" in codon_sequence
            protected = has_cg
            codon = Codon(
                position=i,
                sequence=codon_sequence,
                next_codon=next_codon,
                has_cg=has_cg,
                protected=protected,
                cg_split=cg_split,
            )
            codons.append(codon)

        for codon in codons:
            if codon.has_cg:
                self.original_cg_positions.append(codon.get_cg_position(codon.sequence))
        self.original_cg_positions.sort()
        return codons

    def enforce_packaging_signal(self, codon):
        """
        Ensures a given codon mutation respects the packaging signal constraints.

        There may be more than one mutable position in a codon, so this method checks all of them.

        Parameters:
            codon (Codon): The codon to check for potential mutation.

        Returns:
            bool: True if the mutation respects the packaging signal constraints, False otherwise.
        """

        mutation_positions = codon.get_mutation_positions()
        comparison_results = []

        # Edge case
        if self.packaging_signal_length_end > 1:
            end_length = self.packaging_signal_length_end + 1
        else:
            end_length = self.packaging_signal_length_end

        # True if can mutate, false if can't
        for mutation_position in mutation_positions:
            result = (
                (self.packaging_signal_length_beginning - 1)
                < mutation_position
                < (self.sequence_length - end_length)
            )
            comparison_results.append(result)
        return all(comparison_results)

    def load_new_sequence(self):
        """
        Construct the new nucleotide sequence after mutations have been applied.

        This method updates the new_sequence attribute of the Gene object with the mutated sequence.
        """

        self.new_sequence = ""
        for codon in self.current_codons:
            if codon.mutated:
                self.new_sequence = self.new_sequence + codon.potential_new_seq
            else:
                self.new_sequence = self.new_sequence + codon.sequence

    def determine_changeable_CpG(self):
        """
        Identify codons that can be synonymously mutated to increase CpG sites while respecting gap constraints.

        This method updates the list of mutable positions that can potentially increase CpG sites.
        """

        for codon in self.current_codons:
            if (
                not codon.has_cg
                and not codon.protected
                and (
                    codon.sequence in DIC
                    or (codon.sequence in DIC_for_split and codon.next_codon[0] == "G")
                )
            ):
                if codon.sequence in DIC:
                    codon.potential_new_seq = DIC[codon.sequence]
                else:
                    codon.potential_new_seq = DIC_for_split[codon.sequence]

                cg_position = codon.get_cg_position(codon.potential_new_seq)
                if cg_position is not None:
                    search_start = max(cg_position - self.minimum_CpG_gap_extra, 0)
                    search_end = cg_position + self.minimum_CpG_gap_extra + 2
                    if "CG" not in self.original_sequence[search_start:search_end]:
                        codon.mutable = True
                        self.mutable_positions.append(cg_position)
        self.mutable_positions.sort()

    def translate_sequence(self, sequence):
        """
        Translate a nucleotide sequence into its corresponding amino acid sequence.

        Parameters:
            sequence (str): The nucleotide sequence to translate.

        Returns:
            str: The translated amino acid sequence.
        """

        protein = []
        start = 0
        while start + 2 < len(sequence):
            codon = sequence[start : start + 3]
            protein.append(dna_to_pro[codon][0])
            start += 3
        return "".join(protein)

    def first_difference(self, str1, str2):
        """
        Identify the first difference between two amino acid sequences. Used for debugging purposes.

        Parameters:
            str1 (str): The first amino acid sequence.
            str2 (str): The second amino acid sequence.

        Returns:
            str: The first differing amino acids between the two sequences.
        """

        for a, b in zip(str1, str2):
            if a != b:
                return a + b

    def check_synonymity(self):
        """
        Ensure that mutations applied do not alter the encoded amino acid sequence.

        This method checks for synonymity between the original and new sequences, raising an exception if a non-synonymous mutation is found.
        """

        translated1 = self.translate_sequence(self.original_sequence)
        translated2 = self.translate_sequence(self.new_sequence)

        if translated1 != translated2:
            raise Exception(
                "Code error: Non-synonymous mutations were introduced!\n"
                + self.first_difference(translated1, translated2)
            )

    def check_minimum_gaps(self):
        """
        Validate that the minimum gap constraint between CpG sites is respected.

        This method checks if the mutated sequence maintains the required minimum gaps between CpG sites, raising an exception if the constraint is violated.
        """

        for i in range(1, len(self.current_cg_positions)):
            diff = abs(
                self.current_cg_positions[i]
                - self.current_cg_positions[i - 1]
                - 2  # Subtract 2 because we're looking at the gap
            )
            if diff < self.minimum_CpG_gap and (
                (self.current_cg_positions[i] or self.current_cg_positions[i - 1])
                not in self.original_cg_positions
            ):
                raise Exception(
                    f"Code error: The minimum gap between consecutive CpGs was not "
                    f"preserved! Gap: {diff} at {self.current_cg_positions[i]}"
                )

    def check_packaging_signals(self):
        """
        Validate that mutations do not alter the designated packaging signal regions.

        This method checks for changes in the packaging signal regions due to mutations, raising an exception if such changes are detected.
        """

        if self.packaging_signal_length_end > 0:
            protected_substring1 = (
                self.original_sequence[: self.packaging_signal_length_beginning]
                + self.original_sequence[-self.packaging_signal_length_end :]
            )
            protected_substring2 = (
                self.new_sequence[: self.packaging_signal_length_beginning]
                + self.new_sequence[-self.packaging_signal_length_end :]
            )
        else:
            protected_substring1 = self.original_sequence[
                : self.packaging_signal_length_beginning
            ]
            protected_substring2 = self.new_sequence[
                : self.packaging_signal_length_beginning
            ]

        if protected_substring1 != protected_substring2:
            print()
            print(protected_substring1)
            print(protected_substring2)
            print()
            raise Exception(
                "Code error: Protected packaging signal nucleotides have been changed!"
            )

    def mutate_CpG(self):
        """
        Apply mutations to increase CpG sites while ensuring synonymity and respecting gap and packaging constraints.

        This method applies synonymous mutations to the sequence to increase the number of CpG sites, following the defined constraints.
        """

        for codon in self.current_codons:
            if codon.mutable and all(
                position
                not in range(
                    codon.get_cg_position(codon.potential_new_seq)
                    - self.minimum_CpG_gap_extra,
                    codon.get_cg_position(codon.potential_new_seq)
                    + self.minimum_CpG_gap_extra,
                )
                for position in self.current_cg_positions
            ):
                if self.enforce_packaging_signal(codon):
                    codon.mutated = True
                    self.current_cg_positions.append(
                        codon.get_cg_position(codon.potential_new_seq)
                    )
                    self.current_cg_positions.sort()
        self.load_new_sequence()
        self.check_synonymity()
        self.check_packaging_signals()
        self.check_minimum_gaps()

    def mutate_A_rich(self):
        """
        Mutate the remaining sequence synonymously to increase the prevalence of adenine nucleotides.

        This method targets specific codons for mutation to make the sequence richer in adenine nucleotides while respecting packaging and CpG constraints.
        """

        for codon in self.current_codons:
            if (
                not codon.protected
                and not codon.mutated
                and not codon.has_cg
                and (codon.sequence in DIC_for_A_rich)
            ):
                codon.potential_new_seq = DIC_for_A_rich[codon.sequence]
                if self.enforce_packaging_signal(codon):
                    codon.mutated = True
            # Added this to create A-richer arginine codons, despite protection
            elif (
                codon.sequence == "CGT"
                or codon.sequence == "CGC"
                or codon.sequence == "CGG"
            ):
                codon.potential_new_seq = DIC_for_A_rich[codon.sequence]
                if self.enforce_packaging_signal(codon):
                    codon.mutated = True
        self.load_new_sequence()
        self.check_synonymity()
        self.check_packaging_signals()
        self.check_minimum_gaps()

    def find_desired_gap(self, desired_gap):
        """
        Find the closest possible gap setting to achieve a desired average gap between CpG sites.

        This method iteratively adjusts the minimum CpG gap to approach the desired average gap between CpG sites as closely as possible.

        Parameters:
            desired_gap (int): The target average gap between CpG sites in the sequence.

        Returns:
            int: The closest possible gap setting to the desired average gap.
        """

        def create_sequence(min_gap):
            self.minimum_CpG_gap = min_gap
            self.current_cg_positions = self.original_cg_positions.copy()
            self.current_codons = copy.deepcopy(self.original_codons)
            self.mutable_positions = []
            self.determine_changeable_CpG()
            self.mutate_CpG()

        lower_bound = 0
        upper_bound = self.sequence_length
        closest_gap = lower_bound
        smallest_diff = float("inf")

        while lower_bound <= upper_bound:
            mid = (lower_bound + upper_bound) // 2
            create_sequence(mid)
            avg_gap = self.calculate_average_gap()

            if avg_gap == desired_gap:
                return mid
            elif avg_gap < desired_gap:
                lower_bound = mid + 1
            else:
                upper_bound = mid - 1

            if abs(avg_gap - desired_gap) < smallest_diff:
                smallest_diff = abs(avg_gap - desired_gap)
                closest_gap = mid
        return closest_gap

    def count_CpGs(self, sequence):
        """
        Count the number of CpG sites within the gene sequence.

        A CpG site is identified as a cytosine (C) followed immediately by a guanine (G) in the 5' to 3' direction of the DNA sequence.

        Returns:
            int: The total number of CpG sites found in the gene sequence.
        """

        return sequence.count("CG")

    def calculate_average_gap(self, mode="current", st_dev=False):
        """
        Calculate the average length of gaps between CpG sites in the gene sequence.

        This method identifies all CpG sites within the sequence and calculates the distances between consecutive CpG sites. The average gap length is computed as the mean of these distances.

        Returns:
            float: The average length of gaps between CpG sites in the gene sequence. If there are fewer than two CpG sites, the method returns None, indicating that an average gap cannot be calculated.
        """

        if mode == "original":
            positions = self.original_cg_positions
        elif mode == "current":
            positions = self.current_cg_positions

        gaps = []
        for i in range(len(positions) - 1):
            gap = positions[i + 1] - (positions[i] + 1)
            gaps.append(gap)

        if not gaps:
            return
        else:
            average = round((sum(gaps) / len(gaps)), 3)

        if st_dev is False:
            return average
        if st_dev is True and average is not None:
            variances = [(gap - average) ** 2 for gap in gaps]
            variance_average = sum(variances) / len(variances)
            std_dev = round(math.sqrt(variance_average), 3)
            return std_dev

    def calculate_CpG_abundance_change(self):
        """
        Calculate the percentage change in CpG abundance between the original gene sequence and the mutated gene sequence.

        Returns:
            float: The percentage change in CpG abundance between the new gene sequence and the old sequence. A positive value indicates an increase in CpG abundance from the current sequence to the other, while a negative value indicates a decrease.
        """

        original_CpG_abundance = self.count_CpGs(self.original_sequence)
        final_CpG_abundance = self.count_CpGs(self.new_sequence)

        if original_CpG_abundance == 0:
            if final_CpG_abundance > 0:
                # If original is 0 but final is not, we can consider this as 100%
                # increase
                change = 100.0
            else:
                change = 0.0  # If both original and final are 0, there's no change
        else:
            change = round(
                (
                    (final_CpG_abundance - original_CpG_abundance)
                    / original_CpG_abundance
                )
                * 100,
                1,
            )

        return change

    def calculate_A_abundance_change(self):
        """
        Calculates the percentage change in the abundance of adenine between the original and new gene sequences.

        Returns:
            float: The percentage change in adenine (A) abundance, rounded to one decimal place. A positive value indicates an increase in 'A' abundance, while a negative value indicates a decrease. If the original abundance is 0, returns 0.0.
        """

        original_A_abundance = self.original_sequence.count("A") / len(
            self.original_sequence
        )
        final_A_abundance = self.new_sequence.count("A") / len(self.new_sequence)

        if original_A_abundance == 0:
            return 0.0
        abundance_change = round(
            (((final_A_abundance - original_A_abundance) / original_A_abundance) * 100),
            1,
        )
        return abundance_change


coloring_mode_btn = pnw.RadioButtonGroup(
    name="Coloring mode",
    options={"Highlight CpG's": 1, "Color all nucleotides": 2},
    button_style="solid",
    button_type="light",
)
coloring_mode_btn.visible = False
current_coloring_mode = 1


def view_alignment():
    """
    Generates a visual alignment of sequences using Bokeh for interactive visualization.

    This function extracts the sequences and their IDs from the `sequences` object and prepares them for visualization.
    It employs a nested function, `get_colors`, to assign colors to each base in the sequence depending on the selected coloring mode, highlighting specific features like CG dinucleotides or individual base types.

    The function then constructs two Bokeh plots: one for an overview of the entire sequence alignment (without text) and another detailed view displaying the sequence text with the ability to scroll along the x-axis. The views use color coding for bases and allows for interactive exploration of the sequence alignment.

    Coloring modes can be switched interactively, updating the visualization accordingly. This interactive feature relies on the `pn.depends` decorator to watch for changes in the coloring mode and update the plot data source with new colors.

    Returns:
        Panel object: A Bokeh gridplot object that combines the overview and detailed sequence views, configured for interactive web display.
    """

    global sequences

    if not sequences:  # if list is empty
        return pn.pane.Str("No data")

    # make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (sequences)]
    ids = [rec.id for rec in sequences]
    text = [i for s in seqs for i in s]
    text_str = "".join(text)

    def get_colors(mode):
        # Make colors for bases in sequence. Color 'C' and 'G' green if they form 'CG'
        cache_key = (mode, text_str)
        if cache_key in pn.state.cache:
            return pn.state.cache[cache_key]

        if mode == 1:
            clrs = {
                "A": "white",
                "T": "white",
                "G": "white",
                "C": "white",
                "-": "white",
            }
            colors = [clrs[i] for i in text]
            # Re-color 'C' and 'G' in 'CG' to green
            for i in range(len(colors) - 1):
                if text[i] == "C" and text[i + 1] == "G":
                    colors[i] = "green"  # 'C' in 'CG' is green
                    colors[i + 1] = "green"  # 'G' in 'CG' is green
        elif mode == 2:
            clrs = {"A": "red", "T": "green", "G": "orange", "C": "blue", "-": "white"}
            colors = [clrs[i] for i in text]

        pn.state.cache[cache_key] = colors
        return colors

    colors = get_colors(mode=current_coloring_mode)
    N = len(seqs[0])
    S = len(seqs)

    x = np.arange(1, N + 1)
    y = np.arange(S - 1, -1, -1)
    # creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    # flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    # use recty for rect coords with an offset
    recty = gy + 0.5
    # now we can create the ColumnDataSource with all the arrays

    source = ColumnDataSource(
        data=dict(x=gx, y=gy, recty=recty, text=text, colors=colors)
    )
    plot_height = len(seqs) * 15 + 50
    x_range = Range1d(0, N + 1, bounds="auto", min_interval=50)
    if N > 100:
        viewlen = 100
    else:
        viewlen = N
    # view_range is for the close up view
    view_range = Range1d(0, viewlen, bounds=(0, N + 1))
    tools = "xpan, xwheel_zoom, reset"

    # entire sequence view (no text, with zoom)
    p_head = figure(
        title=None,
        sizing_mode="scale_width",
        height=50,
        x_range=x_range,
        y_range=(0, S),
        tools=tools,
        min_border=0,
        active_scroll="xwheel_zoom",
    )
    rects_head = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        line_width=0,
        fill_alpha=0.6,
        fill_color="colors",
    )
    p_head.add_glyph(source, rects_head)
    p_head.yaxis.visible = False
    p_head.grid.visible = False

    # sequence text view with ability to scroll along x axis
    p_detail = figure(
        title=None,
        sizing_mode="scale_width",
        height=plot_height,
        x_range=view_range,
        y_range=ids[::-1],
        tools="xpan, reset, xwheel_pan",
        min_border=0,
        active_scroll="xwheel_pan",
        lod_threshold=1000000000,  # Prevent bokeh bug when displaying >2000 nucleotides
    )
    glyph = Text(
        x="x",
        y="y",
        text="text",
        text_align="center",
        text_color="black",
        text_font_size="9pt",
    )
    rects_detail = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        fill_color="colors",
        line_width=0,
        fill_alpha=0.4,
    )
    p_detail.add_glyph(source, glyph)
    p_detail.add_glyph(source, rects_detail)

    p_detail.grid.visible = False
    p_detail.xaxis.major_label_text_font_style = "bold"
    p_detail.yaxis.minor_tick_line_width = 0
    p_detail.yaxis.major_tick_line_width = 0

    ticker = SingleIntervalTicker(interval=15, num_minor_ticks=5)
    p_detail.xaxis.ticker = ticker

    p = gridplot([[p_head], [p_detail]], toolbar_location="below")

    @pn.depends(coloring_mode_btn.param.value, watch=True)
    def update_color(change_coloring_mode):
        global current_coloring_mode

        current_coloring_mode = change_coloring_mode
        colors = get_colors(mode=change_coloring_mode)
        source.data.update(colors=colors)

    return p


w1 = pnw.RadioBoxGroup(
    name="Gap method",
    options={"Define minimum gap": 1, "Set target average gap (optimize min gap)": 2},
)
w2 = pnw.IntInput(start=0, placeholder="Your value", sizing_mode="stretch_width")

w3_start = pnw.IntInput(
    start=0, placeholder="Initial length", sizing_mode="stretch_width"
)
w3_end = pnw.IntInput(
    start=0, placeholder="Terminal length", sizing_mode="stretch_width"
)
# Define labels
label_1 = pn.pane.HTML("Initial", styles={"width": "100%", "text-align": "center"})
label_2 = pn.pane.HTML("Final", styles={"width": "100%", "text-align": "center"})

# Define GridBox layout
w3_packaging_signals = pn.GridBox(
    label_1, label_2, w3_start, w3_end, ncols=2, sizing_mode="stretch_width"
)

new_seq_id = pnw.TextInput(
    placeholder="Enter an identifier", sizing_mode="stretch_width"
)
mutate_btn = pnw.Button(name="Mutate", button_type="primary")
A_rich_toggle = pnw.RadioButtonGroup(
    name="A-rich",
    options={"Yes": 1, "No": 2},
    value=2,
    button_style="solid",
    button_type="light",
)

modifiers = pn.FlexBox(
    "# Mutation settings",
    "## CpG gap options",
    w1,
    w2,
    "## Protect terminal nucleotides",
    w3_packaging_signals,
    "## Increase A-richness in remaining codons",
    A_rich_toggle,
    "## Give a unique sequence name",
    new_seq_id,
    pn.FlexBox(mutate_btn, styles={"margin-top": "10px"}, justify_content="center"),
    flex_direction="column",
    styles={"flex-wrap": "wrap"},
)
modifiers.visible = False

# 'How to cite' section will be added here after publication
footer = pn.pane.HTML(
    """

    """,
    styles={
        "width": "100%",
        "text-align": "center",
        "margin-top": "20px",
        "margin-bottom": "20px",
    },
)

template.sidebar.extend([modifiers, footer])

bokeh_pane = pn.pane.Bokeh(sizing_mode="stretch_width", margin=10)

sequences = []


def download_alignment():
    """
    Converts the in-memory sequence alignment into a downloadable FASTA format.

    This function takes the global `sequences` variable, which is expected to be a list of Bio.SeqRecord objects
    (from Biopython), and writes these sequences to a FASTA formatted string. This string is then encoded to bytes,
    which can be used for creating a file-like object suitable for downloading.

    Returns:
        BytesIO: A file-like object containing the sequences in FASTA format, encoded as bytes.
    """

    output = StringIO()
    SeqIO.write(sequences, output, "fasta")
    fasta_str = output.getvalue()
    fasta_bytes = fasta_str.encode()

    return BytesIO(fasta_bytes)


download_btn = pnw.FileDownload(
    callback=download_alignment,
    filename="alignment.fasta",
    button_style="solid",
    button_type="primary",
    height=32,
    label="Download alignment",
)
download_btn.visible = False
successful_load_dummy = pnw.Checkbox(value=False, visible=False)


def load_sample(event):
    """
    Load the sample sequence into the application and update the UI.

    This function is triggered when the user clicks the "Example sequence" button. It loads the sample sequence into the
    application, resets the global state, and updates the UI to reflect the new sequence.

    Parameters:
        event (bokeh.events.ButtonClick): The event object triggered by the button click.
    """

    global sequences, current_coloring_mode, df, download_btn

    sequences = []
    # Reset the dataframe
    df = pd.DataFrame(
        columns=[
            "Name",
            "Average CpG gap",
            "Mutation settings",
            "CpG count",
            "CpG change",
            "A change",
        ]
    )
    df.index.name = "Index"
    df.index = df.index + 1

    sample_string = (
        "AGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACACATAATCCACCTATCC"
        "CAGTAGGAGAAATCTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCAT"
        "TCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGATTCTATAAAACTCTAAGAGCCGAG"
        "CAAGCTTCACAAGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTA"
        "TTTTAAAAGCATTGGGACCAGGAGCGACACTAGAAGAAATGATGACAGCATGTCAGGGAGTGGGGGGACCCGGCCA"
        "TAAAGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATCCAGCTACCATAATGATACAGAAAGGCAATTTT"
        "AGGAACCAAAGAAAGACTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACATAGCC"
    )
    sample_record = SeqRecord(
        Seq(sample_string),
        id="Gag",
        name="Gag",
        description="Sample Gag sequence of HIV-1",
    )

    def encode_fasta_as_bytes(seq_record):
        fasta_stream = StringIO()
        SeqIO.write(seq_record, fasta_stream, "fasta")
        fasta_content = fasta_stream.getvalue()
        return fasta_content.encode("utf-8")

    file_input.value = encode_fasta_as_bytes(sample_record)


load_sample_btn = pnw.Button(
    name="Example sequence",
    button_type="default",
    styles={"margin-top": "1em", "margin-bottom": "1em"},
)
load_sample_btn.on_click(load_sample)

file_input = pnw.FileInput(accept=".fasta,.fa,.fas,.aln")


# Define function to load FASTA file and display sequences
def load_fasta(event):
    """
    Load a FASTA file and display the sequences in the application.

    This function is triggered when a new FASTA file is uploaded using the file input widget. It reads the file
    contents, parses the sequences, and displays the alignment in the application. It also updates the global state
    and resets the mutation settings.

    Parameters:
        event (bokeh.events.ButtonClick): The event object triggered by the file input widget.
    """

    global sequences, current_coloring_mode, df
    global top, download_btn, file_input_watcher, file_input

    sequences = []
    # Reset the dataframe
    df = pd.DataFrame(
        columns=[
            "Name",
            "Average CpG gap",
            "Mutation settings",
            "CpG count",
            "CpG change",
            "A change",
        ]
    )
    df.index.name = "Index"
    df.index = df.index + 1

    fasta_bytes = event.new
    # Get the loaded file bytes
    if fasta_bytes:
        # Decode the bytes to a string
        fasta_string = fasta_bytes.decode("utf-8")
        # Create a StringIO object from the string
        fasta_io = StringIO(fasta_string)
        # Parse the sequences
        sequences = list(SeqIO.parse(fasta_io, "fasta"))

        sequence_lengths = [len(sequence) for sequence in sequences]
        if len(sequences) > 1 and len(set(sequence_lengths)) > 1:
            top_message.visible = True
            notifications.warning(
                "The sequences in your input file must be of equal length!", duration=0
            )
            return

        for sequence in sequences:
            sequence.seq = sequence.seq.upper()
            # Turn U's into T's for convenience
            sequence.seq = sequence.seq.back_transcribe()

        # Add data to the DataFrame for each sequence
        for sequence in sequences:
            try:
                gene = Gene(sequence.seq, gap_method=None)
            except Exception as e:
                notifications.warning(str(e), duration=0)
                top_message.visible = True
                return
            else:
                update_table(
                    average_CpG_gap=gene.original_average_gap,
                    sequence_id=sequence.id,
                    CpG_count=gene.count_CpGs(gene.original_sequence),
                )

        # Display the sequences
        p = view_alignment()
        bokeh_pane.object = p

        modifiers.visible = True
        top_message.visible = False

        if coloring_mode_btn.visible is False:
            coloring_mode_btn.visible = True

        download_btn.filename = f"{sequences[0].id}-recoding-aln.fasta"
        download_btn.visible = len(sequences) > 1

        successful_load_dummy.value = not successful_load_dummy.value

        file_input.param.unwatch(file_input_watcher)
        new_file_input = pnw.FileInput(accept=".fasta,.fa,.fas,.aln")
        new_file_input_watcher = new_file_input.param.watch(load_fasta, "value")
        top[1] = new_file_input
        file_input = new_file_input
        file_input_watcher = new_file_input_watcher


file_input_watcher = file_input.param.watch(load_fasta, "value")


def mutate(
    event,
    gap_method,
    option_value,
    protection_start_value,
    protection_end_value,
    identifier,
    A_rich,
):
    """
    Apply synonymous mutations to the sequence according to the user's settings.

    This function applies synonymous mutations to the sequence according to the user's settings, including the gap method,
    gap value, protection settings, and A-richness. It then updates the sequence alignment display and the mutation
    settings table.

    Parameters:
        event (bokeh.events.ButtonClick): The event object triggered by the mutation button.
        gap_method (int): The method for specifying the minimum gap between CpG sites. If 1, the user specifies the
            minimum gap value. If 2, the user specifies the desired average gap value, and the algorithm optimizes the
            minimum gap value.
        option_value (int): The minimum or desired average gap value between CpG sites.
        protection_start_value (int): The number of initial nucleotides to protect from mutation.
        protection_end_value (int): The number of terminal nucleotides to protect from mutation.
        identifier (str): The identifier for the mutated sequence.
        A_rich (bool): Whether to increase A-richness in the remaining codons.
    """

    global sequences, download_btn

    original_id = sequences[0].id
    description = ""
    mutation_settings = ""

    if gap_method == 1:
        try:
            gene = Gene(
                sequences[0].seq,
                gap_method,
                minimum_CpG_gap=option_value,
                packaging_signal_length_beginning=protection_start_value,
                packaging_signal_length_end=protection_end_value,
            )
        except Exception as e:
            notifications.error(str(e), duration=0)
            return
        description += (
            f"Synonymously recoded {original_id} sequence "
            f"to increase CpG's with a minimum gap of {option_value}"
        )
        mutation_settings += f"Minimum CpG gap set as {option_value}. "
    elif gap_method == 2:
        try:
            gene = Gene(
                sequences[0].seq,
                gap_method,
                minimum_CpG_gap=None,
                desired_CpG_gap=option_value,
                packaging_signal_length_beginning=protection_start_value,
                packaging_signal_length_end=protection_end_value,
            )
        except Exception as e:
            notifications.error(str(e), duration=0)
            return
        resultant_avg_gap = round(gene.calculate_average_gap())
        description += (
            f"Synonymously recoded {original_id} sequence "
            f"to increase CpG's at an average gap of {resultant_avg_gap}"
        )
        mutation_settings += (
            f"Desired average CpG gap set as {option_value}.\n"
            f"Algorithm-optimized minimal gap value is {gene.minimum_CpG_gap}. "
        )
    if A_rich:
        try:
            gene.mutate_A_rich()
        except Exception as e:
            notifications.error(str(e), duration=0)
            return
        description += ". Mutated remaining codons to A-rich versions, if possible."
        mutation_settings += "\nA-enriched remaining codons. "

    sequences.append(
        SeqRecord(
            gene.new_sequence, id=identifier, name=identifier, description=description
        )
    )

    if protection_start_value > 0:
        mutation_settings += (
            f"\nProtected {protection_start_value} initial nucleotides. "
            if protection_end_value == 0
            else (
                f"\nProtected {protection_start_value} initial and "
                f"{protection_end_value} terminal nucleotides. "
            )
        )
    elif protection_end_value > 0:
        mutation_settings += (
            f"\nProtected {protection_end_value} terminal nucleotides. "
        )

    bokeh_pane.object = view_alignment()
    download_btn.visible = len(sequences) > 1

    CpG_change = (
        f"+{gene.calculate_CpG_abundance_change()}"
        if gene.calculate_CpG_abundance_change() > 0
        else gene.calculate_CpG_abundance_change()
    )
    A_change = (
        f"+{gene.calculate_A_abundance_change()}"
        if gene.calculate_A_abundance_change() > 0
        else gene.calculate_A_abundance_change()
    )

    update_table(
        average_CpG_gap=gene.calculate_average_gap(),
        sequence_id=identifier,
        CpG_count=gene.count_CpGs(gene.new_sequence),
        mutation_settings=mutation_settings.strip(),
        CpG_abundance_change=f"{CpG_change}%",
        A_abundance_change=f"{A_change}%",
    )


def on_mutate_btn_click(event):
    """
    Handle the click event for the mutation button.

    This function is triggered when the user clicks the mutation button. It reads the user's settings for the mutation
    and applies the synonymous mutations to the sequence accordingly.

    Parameters:
        event (bokeh.events.ButtonClick): The event object triggered by the mutation button.
    """

    global sequences

    gap_method = w1.value
    option_value = w2.value

    if w2.value is None:
        notifications.warning("Please enter the mutation value.")
        return

    if w3_start.value is None:
        w3_start.value = 0
    if w3_end.value is None:
        w3_end.value = 0
    protection_start_value = w3_start.value
    protection_end_value = w3_end.value
    identifier = new_seq_id.value_input.strip()

    if A_rich_toggle.value == 1:
        A_rich = True
    elif A_rich_toggle.value == 2:
        A_rich = False

    if identifier is None or identifier == "":
        notifications.warning("Please enter at least one character.")
    else:
        found = False
        for seq in sequences:
            if seq.id == identifier:
                found = True
                notifications.warning("A sequence with this name already exists.")
                return

        if not found:
            mutate(
                event=event,
                gap_method=gap_method,
                option_value=option_value,
                protection_start_value=protection_start_value,
                protection_end_value=protection_end_value,
                identifier=identifier,
                A_rich=A_rich,
            )


top = pn.FlexBox(
    "## Load a new FASTA File",
    file_input,
    load_sample_btn,
    successful_load_dummy,
    justify_content="space-evenly",
)
top_message = pn.pane.Markdown(
    """
    ### This app applies synonymous mutations to add CpG dinucleotides into DNA, \
        according to user-specified constraints.
    It can also make the remaining sequence more A-rich, \
        while maintaining the amino acid sequence.

    Please load a sequence starting in-frame, of a nucleotide length divisible by 3. \
        It must only contain codons present in the \
        <a href="https://www.hgmd.cf.ac.uk/docs/cd_amino.html">codon table</a>.

    If loading a FASTA file with more than one sequence, ensure that \
        sequence lengths are equal.

    In such case, the sequence at the top of the file will be used by the \
        recoding algorithm.
    """,
    max_width=600,
)
accessory = pn.FlexBox(coloring_mode_btn, download_btn, justify_content="space-between")

visualization = pn.FlexBox(
    top,
    top_message,
    bokeh_pane,
    accessory,
    flex_direction="column",
    align_items="center",
    sizing_mode="stretch_width",
)


def update_table(
    average_CpG_gap,
    sequence_id,
    CpG_count,
    mutation_settings=None,
    CpG_abundance_change=None,
    A_abundance_change=None,
):
    """
    Update the mutation settings table with new data.

    This function updates the mutation settings table with new data, including the sequence identifier, average CpG gap,
    mutation settings, CpG count, and changes in CpG and A abundance.

    Parameters:
        average_CpG_gap (float): The average gap between CpG sites in the sequence.
        sequence_id (str): The identifier for the sequence.
        mutation_settings (str, optional): The settings used for the synonymous mutations. Defaults to None.
        CpG_count (int): The number of CpG sites in the sequence.
        CpG_abundance_change (float, optional): The percentage change in CpG abundance between the original and new
            sequences. Defaults to None.
        A_abundance_change (float, optional): The percentage change in A abundance between the original and new sequences.
            Defaults to None.
    """

    global df

    data = []

    data = pd.DataFrame(
        [
            {
                "Name": sequence_id,
                "Average CpG gap": average_CpG_gap,
                "Mutation settings": mutation_settings,
                "CpG count": CpG_count,
                "CpG change": CpG_abundance_change,
                "A change": A_abundance_change,
            }
        ]
    )

    # Concatenate new data to existing DataFrame
    df = pd.concat([df, data], ignore_index=True)
    df.index.name = "Index"
    df.index = df.index + 1
    df.fillna("NA", inplace=True)
    tabulator.value = df


df = pd.DataFrame(
    columns=[
        "Name",
        "Average CpG gap",
        "Mutation settings",
        "CpG count",
        "CpG change",
        "A change",
    ]
)
df.index.name = "Index"
df.index = df.index + 1
formatter = {"Mutation settings": {"type": "textarea"}}
tabulator = pnw.Tabulator(df, formatters=formatter, layout="fit_data", disabled=True)

global tableDownloadAdded
tableDownloadAdded = False

table = pn.FlexBox(
    tabulator,
    flex_direction="column",
    align_content="center",
    sizing_mode="stretch_width",
)


template.main.extend(
    [
        visualization,
        table,
    ]
)

successful_load_dummy.jscallback(
    value="""
        if (document.getElementById('sidebar-button') &&
        document.getElementById('sidebar').classList.contains('hidden')) {
            document.getElementById('sidebar-button').click();
        }
        """
)

mutate_btn.param.watch(on_mutate_btn_click, "clicks")

template.servable()
