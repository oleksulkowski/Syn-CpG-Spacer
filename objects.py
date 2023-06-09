import sys, copy, math

# Taken from github.com/Mmark94/protein_recoding. For STOP codons, X is used as symbol.
dna_to_pro = {'ATG': 'M', 'GCG': 'A', 'TCA': 'S', 'GAA': 'E', 'GGG': 'G', 'GGT': 'G', 'AAA': 'K', 'GAG': 'E', 'AAT': 'N', 'CTA': 'L',
                  'CAT': 'H', 'TCG': 'S', 'TAG': 'X', 'GTG': 'V', 'TAT': 'Y', 'CCT': 'P', 'ACT': 'T', 'TCC': 'S', 'CAG': 'Q', 'CCA': 'P',
                  'TAA': 'X', 'AGA': 'R', 'ACG': 'T', 'CAA': 'Q', 'TGT': 'C', 'GCT': 'A', 'TTC': 'F', 'AGT': 'S', 'ATA': 'I', 'TTA': 'L',
                  'CCG': 'P', 'ATC': 'I', 'TTT': 'F', 'CGT': 'R', 'TGA': 'X', 'GTA': 'V', 'TCT': 'S', 'CAC': 'H', 'GTT': 'V', 'GAT': 'D',
                  'CGA': 'R', 'GGA': 'G', 'GTC': 'V', 'GGC': 'G', 'TGC': 'C', 'CTG': 'L', 'CTC': 'L', 'CGC': 'R', 'CGG': 'R', 'AAC': 'N',
                  'GCC': 'A', 'ATT': 'I', 'AGG': 'R', 'GAC': 'D', 'ACC': 'T', 'AGC': 'S', 'TAC': 'Y', 'ACA': 'T', 'AAG': 'K', 'GCA': 'A',
                  'TTG': 'L', 'CCC': 'P', 'CTT': 'L', 'TGG': 'W'}

DIC = {
    "TCT": "TCG", "TCC": "TCG", "TCA": "TCG", "AGT": "TCG", "AGC": "TCG", # Serine
    "CCT": "CCG", "CCC": "CCG", "CCA": "CCG", # Proline
    "ACT": "ACG", "ACC": "ACG", "ACA": "ACG", # Threonine
    "GCT": "GCG", "GCC": "GCG", "GCA": "GCG", # Alanine
    "AGA": "CGA", "AGG": "CGA", # Arginine - don't need to add the rest as they already have a CpG
    }

DIC_for_split = {
    "TTG": "CTC", "TTA": "CTC", "CTT": "CTC", "CTA": "CTC", "CTG": "CTC", # Leucine
    "GTT": "GTC", "GTA": "GTC", "GTG": "GTC", # Valine
    "GAT": "GAC", # Aspartic acid
    "GGT": "GGC", "GGA": "GGC", "GGG": "GGC", # Glycine
    "AGT": "AGC", "TCA": "TCC", "TCT": "TCC", # Serine
    "CGG": "CGC", "CGA": "CGC", "CGT":"CGC", "AGG": "CGC", "AGA": "CGC", # Arginine
    "AAT": "AAC", # Asparagine
    "ATA": "ATC", "ATT": "ATC", # Isoleucine
    "ACG": "ACC", "ACA": "ACC", "ACT": "ACC", # Threonine
    "TGT": "TGC", # Cysteine
    "TAT": "TAC", # Tyrosine
    "TTT": "TTC", # Phenylalanine
    "CAT": "CAC", # Histidine
    "CCA": "CCC", "CCT": "CCC" # Proline
}


DIC_for_A_rich = {
    "TTG": "TTA", "CTT": "CTA", "CTC": "CTA", "CTG": "CTA", # Leucine
    "ATT": "ATA", "ATC": "ATA", # Isoleucine
    "GTT": "GTA", "GTC": "GTA", "GTG": "GTA", # Valine
    "TCT": "TCA", "TCC": "TCA", # Serine
    "CCT": "CCA", "CCC": "CCA", # Proline
    "ACT": "ACA", "ACC": "ACA", # Threonine
    "GCT": "GCA", "GCC": "GCA", # Alanine
    "CAG": "CAA", # Glutamine
    "GAG": "GAA", # Glutamic acid
    "CGT": "CGA", "CGC": "CGA", "CGG": "CGA", # Arginine - fixed the code to allow A-enrichment despite CpG, if it is not in the packaging signal
    "AGG": "AGA", # Arginine
    "GGT": "GGA", "GGC": "GGA", "GGG": "GGA" # Glycine
    # STOP codons are not mutated

}

class Codon:
    def __init__(self, position, sequence, next_codon=None, has_cg=False, protected=False, cg_split=False, mutable=None, potential_new_seq=None, mutated=False):
        self.position = position # Position of the codon, not the CG!
        self.sequence = sequence
        self.next_codon = next_codon
        self.has_cg = has_cg
        self.protected = protected
        self.cg_split = cg_split
        self.mutable = mutable
        self.potential_new_seq = potential_new_seq
        self.mutated = mutated

        self.check_codon_viability()

    # Relative to the 'C'
    def get_cg_position(self, sequence):
        if sequence is not None:
            # Check for a CG pair within this codon
            for i in range(3):
                if sequence[i] == 'C' and i < 2 and sequence[i + 1] == 'G':
                    return self.position + i

            # Check for a CG pair split across two codons
            if (
                self.next_codon is not None
                and sequence[2] == 'C'
                and self.next_codon[0] == 'G'
            ):
                return self.position + 2

        # No CG pair was found
        return None

    def get_mutation_positions(self):
        diff_positions = []

        for i in range(3):
            if self.sequence[i] != self.potential_new_seq[i]:
                diff_positions.append(self.position + i)
        return diff_positions

    def translate_codon(self):
        if self.sequence in dna_to_pro:
            aa_symbol = dna_to_pro[self.sequence]
            return aa_symbol
    
    # Checks if the codon is viable. If not, this function will raise an error and exit the program
    def check_codon_viability(self):
        aa_symbol = self.translate_codon()
        if aa_symbol is None:
            sys.exit('The file you loaded contains invalid codons/nucleotides. Please ensure your sequence contains full, viable, in-frame codons')


class Gene:
    def __init__(self, original_sequence, gap_method, packaging_signal_length_beginning = 0, packaging_signal_length_end = 0, minimum_CpG_gap = 12, desired_CpG_gap = None):
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

        self.original_codons = self.analyze_codons() # Used only to store the original codons - no operations are performed on this attribute
        self.current_codons = copy.deepcopy(self.original_codons)
        self.current_cg_positions = self.original_cg_positions.copy()
        self.original_average_gap = self.calculate_average_gap("original")

        self.check_gap_method()

    @property
    def minimum_CpG_gap_extra(self):
        return self.minimum_CpG_gap + 1

    # See if the user defined a minimum gap, or defined an average CpG spacing, towards which the gap will be adjusted
    def check_gap_method(self):
        if self.gap_method == None:
            return
        elif self.gap_method == 1 and self.desired_CpG_gap == None:
            self.determine_changeable_CpG()
            self.mutate_CpG()
        elif self.gap_method == 2 and self.minimum_CpG_gap == None:
            self.closest_gap = self.find_desired_gap(self.desired_CpG_gap)
        elif self.gap_method == 3 and self.minimum_CpG_gap == None and self.desired_CpG_gap == None:
            def prompt_for_integer(message):
                while True:
                    try:
                        user_input = int(input(message))
                        if user_input < 0:
                            print("This is a negative value. Please try again.")
                        else:
                            break
                    except ValueError:
                        print("This is not a valid integer. Please try again.")
                
                return user_input
            
            print("The value will be found between your specified average gap length bounds")
            min_avg_bound = prompt_for_integer('Please enter the minimum average gap value (integer): ')
            max_avg_bound = prompt_for_integer('Now, please enter the maximum value: ')
            self.find_min_cv_gap(min_avg_bound, max_avg_bound)

    # Creates class instances for codons and looks for original CpG's
    def analyze_codons(self):
        codons = []
        for i in range(0, self.sequence_length, 3):
            codon_sequence = self.original_sequence[i:i + 3]
            has_cg = 'CG' in codon_sequence
            cg_split = None
            if i + 3 < len(self.original_sequence):
                next_codon = self.original_sequence[i + 3:i + 6]
                if 'CG' in (codon_sequence[-1] + next_codon[0]):
                    has_cg = True
                    cg_split = True
            else:
                has_cg = 'CG' in codon_sequence
            protected = has_cg
            codon = Codon(position=i, sequence=codon_sequence, next_codon=next_codon, has_cg=has_cg, protected=protected, cg_split=cg_split)
            codons.append(codon)
        
        for codon in codons:
            if codon.has_cg:
                self.original_cg_positions.append(codon.get_cg_position(codon.sequence))
        self.original_cg_positions.sort()
        return codons
    
    def enforce_packaging_signal(self, codon):
        mutation_positions = codon.get_mutation_positions()
        comparison_results = []
        for mutation_position in mutation_positions:
                result = (self.packaging_signal_length_beginning - 1) < mutation_position < (self.sequence_length - self.packaging_signal_length_end)
                comparison_results.append(result)
                return all(comparison_results)

    # Applies mutations into the object to create a new sequence
    def load_new_sequence(self):
        self.new_sequence = ''
        for codon in self.current_codons:
            if codon.mutated:
                self.new_sequence = self.new_sequence + codon.potential_new_seq
            else:
                self.new_sequence = self.new_sequence + codon.sequence

    # Finds which codons can be synonymously mutated to give more CpG's. Such codons must have at least a minimum_CpG_gap nucleotides long gap between another CpG
    def determine_changeable_CpG(self):
        for codon in self.current_codons:
            if not codon.has_cg and not codon.protected and (codon.sequence in DIC or (codon.sequence in DIC_for_split and codon.next_codon[0] == 'G')):
                if codon.sequence in DIC:
                    codon.potential_new_seq = DIC[codon.sequence]
                else:
                    codon.potential_new_seq = DIC_for_split[codon.sequence]
                
                for i in range(self.sequence_length):
                    if i <= self.minimum_CpG_gap and self.original_sequence[:i+(self.minimum_CpG_gap_extra+2)].find('CG') == -1 and codon.get_cg_position(codon.potential_new_seq) == i:
                        codon.mutable = True
                    if i > self.minimum_CpG_gap and self.original_sequence[i-(self.minimum_CpG_gap_extra):i+(self.minimum_CpG_gap_extra+2)].find('CG') == -1 and codon.get_cg_position(codon.potential_new_seq) == i:
                        codon.mutable = True
            if codon.mutable:
                self.mutable_positions.append(codon.get_cg_position(codon.potential_new_seq))
        self.mutable_positions.sort()

    # Applies synonymous mutations in a 5' to 3' direction. Ensures there are sufficient gaps between new CpG's
    def mutate_CpG(self):
        for codon in self.current_codons:

            if codon.mutable and all(
                    position not in range(
                        codon.get_cg_position(codon.potential_new_seq) - self.minimum_CpG_gap_extra,
                        codon.get_cg_position(codon.potential_new_seq) + self.minimum_CpG_gap_extra)
                for position in self.current_cg_positions):
                    if (self.enforce_packaging_signal(codon)):
                        codon.mutated = True
                        self.current_cg_positions.append(codon.get_cg_position(codon.potential_new_seq))
                        self.current_cg_positions.sort()
        self.load_new_sequence()

    # Makes the new sequence richer in A nucleotides in codons that wouldn't interfere with packaging or CpG's
    def mutate_A_rich(self):
        for codon in self.current_codons:
            if not codon.protected and not codon.mutated and not codon.has_cg and (codon.sequence in DIC_for_A_rich):
                codon.potential_new_seq = DIC_for_A_rich[codon.sequence]
                if (self.enforce_packaging_signal(codon)):
                        codon.mutated = True
            # Added this to create A-richer arginine codons, though they are meant to be protected
            elif codon.sequence == "CGT" or codon.sequence == "CGC" or codon.sequence == "CGG":
                codon.potential_new_seq = DIC_for_A_rich[codon.sequence]
                if (self.enforce_packaging_signal(codon)):
                        codon.mutated = True
        self.load_new_sequence()

    def find_desired_gap(self, desired_gap, no_print = False):
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
        smallest_diff = float('inf')

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
        
        if no_print == False:
            print(f"A minimum gap of {closest_gap} was found to result in the closest match to your desired gap")
        return closest_gap

    def find_min_cv_gap(self, min_avg_bound, max_avg_bound):
        def create_sequence(min_gap):
            self.minimum_CpG_gap = min_gap
            self.current_cg_positions = self.original_cg_positions.copy()
            self.current_codons = copy.deepcopy(self.original_codons)
            self.mutable_positions = []
            self.determine_changeable_CpG()
            self.mutate_CpG()

        best_gap = max_avg_bound
        smallest_cv = float('inf')


        for gap in range(min_avg_bound, max_avg_bound + 1):
            create_sequence(gap)

            std_gap = self.calculate_average_gap("current", st_dev=True)
            avg_gap = self.calculate_average_gap()

            # calculate the coefficient of variation
            if avg_gap != 0:  # Avoid division by zero
                cv = std_gap / avg_gap
            else:
                cv = float('inf')

            print(f"{gap} {avg_gap} {cv} {self.count_CpGs(self.new_sequence)}")
            if min_avg_bound <= avg_gap <= max_avg_bound and cv < smallest_cv:
                smallest_cv = cv
                best_gap = gap
        

        print(f"A minimum gap of {best_gap} was found to result in the smallest coefficient of variation of the gap, {cv}.")
        return best_gap


    # Translates a sequence into amino acids
    def translate_sequence(self, sequence):
        protein = []
        start = 0
        while start + 2 < len(sequence):
            codon = sequence[start:start + 3]
            protein.append(dna_to_pro[codon][0])
            start += 3
        return ''.join(protein)

    # If there is a difference between the original amino acid sequence and the new one, returns the first difference in amino acids that it finds
    def first_difference(self, str1, str2):
        for a, b in zip(str1, str2):
            if a != b:
                return a+b

    def check_synonymity(self):
        translated1 = self.translate_sequence(self.original_sequence)
        translated2 = self.translate_sequence(self.new_sequence)

        # Checks if the original amino acid sequence is different to the new one and prints an error if so. If original sequence has a stop codon change translated1 for translated1[:-3]
        if translated1 != translated2:
            sys.exit('ERROR: Non-synonymous mutations were introduced!\n' + self.first_difference(translated1, translated2))
    
    # Check if the packaging signals have been modified
    def check_packaging_signals(self):
        if self.packaging_signal_length_end > 0:
            protected_substring1 = self.original_sequence[:self.packaging_signal_length_beginning] + self.original_sequence[-self.packaging_signal_length_end:]
            protected_substring2 = self.new_sequence[:self.packaging_signal_length_beginning] + self.new_sequence[-self.packaging_signal_length_end:]
        else:
            protected_substring1 = self.original_sequence[:self.packaging_signal_length_beginning]
            protected_substring2 = self.new_sequence[:self.packaging_signal_length_beginning]

        if protected_substring1 != protected_substring2:
            sys.exit('ERROR: Protected packaging signal nucleotides have been changed!')

    # Counts the number of CpG's in a genetic sequence
    def count_CpGs(self, sequence):
        counter = 0
        for i in range(len(sequence) - 1):
            if sequence[i:i+2] == "CG":
                counter += 1
        return counter

    # Calculates the average gap between CpG's in the sequence
    def calculate_average_gap(self, mode = "current", st_dev = False):
        if mode == "original":
            positions = self.original_cg_positions
        elif mode == "current":
            positions = self.current_cg_positions

        gaps = []
        for i in range(len(positions) - 1):
            gap = positions[i + 1] - (positions[i] + 1)
            gaps.append(gap)
        
        if len(gaps) == 0:
            return None
        else:
            average = round((sum(gaps) / len(gaps)), 3)

        if st_dev == False:
            return average
        if st_dev == True and average is not None:
            variances = [(gap - average) ** 2 for gap in gaps]
            variance_average = sum(variances) / len(variances)
            std_dev = round(math.sqrt(variance_average), 3)
            return std_dev

    def calculate_CpG_abundance_change(self):
        original_CpG_abundance = self.count_CpGs(self.original_sequence)
        final_CpG_abundance = self.count_CpGs(self.new_sequence)

        if original_CpG_abundance == 0:
            if final_CpG_abundance > 0:
                change = 100.0  # If original is 0 but final is not, we can consider this as 100% increase
            else:
                change = 0.0  # If both original and final are 0, there's no change
        else:
            change = round(((final_CpG_abundance - original_CpG_abundance) / original_CpG_abundance)*100, 1)

        print(f"CpG abundance changed {change}%")
        return change

    # Calculates the change in A-richness of the sequence
    def calculate_A_abundance_change(self):
        original_A_abundance = self.original_sequence.count('A') / len(self.original_sequence)
        final_A_abundance = self.new_sequence.count('A') / len(self.new_sequence)
        print(f"Initial vs final A proportion: {round(original_A_abundance, 3)} vs {round(final_A_abundance, 3)}")

        if original_A_abundance == 0:
            return 0.0
        abundance_change = round((((final_A_abundance - original_A_abundance) / original_A_abundance)*100), 1)
        print(f"A-nucleotide abundance change: {abundance_change}%")
        return abundance_change

    # Unused by default. Prints the original sequence in the terminal
    def read_codons(self):
        print("Original sequence: ", end="")
        for codon in self.original_codons:
            print(codon.sequence, end="")
        print("")

    # Unused by default. Shows the details of which codons can possibly be mutated. Prints into terminal. Good for testing without saving
    def print_changeable_CpG(self):
        for codon in self.current_codons:
            if codon.mutable:
                print(f"{codon.sequence} | {codon.position} | {codon.potential_new_seq} | {codon.get_cg_position(codon.potential_new_seq)}")
        print("Codon | Codon position | New codon | New CpG position (based on C)")
        print("")