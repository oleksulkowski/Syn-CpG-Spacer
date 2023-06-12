importScripts("https://cdn.jsdelivr.net/pyodide/v0.23.0/full/pyodide.js");

function sendPatch(patch, buffers, msg_id) {
  self.postMessage({
    type: 'patch',
    patch: patch,
    buffers: buffers
  })
}

async function startApplication() {
  console.log("Loading pyodide!");
  self.postMessage({type: 'status', msg: 'Loading pyodide'})
  self.pyodide = await loadPyodide();
  self.pyodide.globals.set("sendPatch", sendPatch);
  console.log("Loaded!");
  await self.pyodide.loadPackage("micropip");
  const env_spec = ['markdown-it-py<3', 'https://cdn.holoviz.org/panel/1.1.0/dist/wheels/bokeh-3.1.1-py3-none-any.whl', 'https://cdn.holoviz.org/panel/1.1.0/dist/wheels/panel-1.1.0-py3-none-any.whl', 'pyodide-http==0.2.1', 'Bio', 'numpy', 'pandas']
  for (const pkg of env_spec) {
    let pkg_name;
    if (pkg.endsWith('.whl')) {
      pkg_name = pkg.split('/').slice(-1)[0].split('-')[0]
    } else {
      pkg_name = pkg
    }
    self.postMessage({type: 'status', msg: `Installing ${pkg_name}`})
    try {
      await self.pyodide.runPythonAsync(`
        import micropip
        await micropip.install('${pkg}');
      `);
    } catch(e) {
      console.log(e)
      self.postMessage({
	type: 'status',
	msg: `Error while installing ${pkg_name}`
      });
    }
  }
  console.log("Packages loaded!");
  self.postMessage({type: 'status', msg: 'Executing code'})
  const code = `
  
import asyncio

from panel.io.pyodide import init_doc, write_doc

init_doc()

import pandas as pd
import numpy as np
from io import StringIO, BytesIO
import sys, copy, math

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

import panel as pn
import panel.widgets as pnw

pn.extension('tabulator', 'notifications', design='fast', notifications=True)
template = pn.template.FastListTemplate(title='SynRecoder', theme_toggle=False, accent='#A01346', collapsed_sidebar=True, main_max_width='80vw')

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
            sys.exit('ERROR: Non-synonymous mutations were introduced!\\n' + self.first_difference(translated1, translated2))
    
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

coloring_mode_btn = pnw.RadioButtonGroup(name='Coloring mode', options={"Highlight CpG's": 1, 'Color all nucleotides': 2}, button_style='solid', button_type='light')
current_coloring_mode = 1
seqs = pnw.LiteralInput(value=[], type=list)
@pn.depends(seqs.param.value, coloring_mode_btn.param.value)
def view_alignment():
    """Bokeh sequence alignment view"""
    global sequences

    if not sequences: # if list is empty
        return pn.pane.Str('No data')

    #make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (sequences)]
    ids = [rec.id for rec in sequences]
    text = [i for s in seqs for i in s]
    colors = get_colors(seqs)    
    N = len(seqs[0])
    S = len(seqs)    
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(S-1,-1,-1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5
    #now we can create the ColumnDataSource with all the arrays
    
    source = ColumnDataSource(data=dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs)*15+50
    x_range = Range1d(0,N+1, bounds='auto', min_interval=50)
    if N>100:
        viewlen=100
    else:
        viewlen=N
    #view_range is for the close up view
    view_range = Range1d(0, viewlen, bounds=(0, N+1))
    tools="xpan, xwheel_zoom, reset"

    #entire sequence view (no text, with zoom)
    p_head = figure(title=None, sizing_mode='scale_width', height=50,
               x_range=x_range, y_range=(0, S), tools=tools,
               min_border=0, active_scroll='xwheel_zoom')
    rects_head = Rect(x="x", y="recty", width=1, height=1, line_width=0, fill_alpha=0.6, fill_color="colors")
    p_head.add_glyph(source, rects_head)
    p_head.yaxis.visible = False
    p_head.grid.visible = False

    #sequence text view with ability to scroll along x axis
    p_detail = figure(title=None, sizing_mode='scale_width', height=plot_height,
                x_range=view_range, y_range=ids[::-1], tools="xpan, reset, xwheel_pan",
                min_border=0, active_scroll='xwheel_pan', lod_threshold=1000000000)          
    glyph = Text(x="x", y="y", text="text", text_align="center", text_color="black",
                text_font_size="9pt")
    rects_detail = Rect(x="x", y="recty", width=1, height=1, fill_color="colors", line_width=0, fill_alpha=0.4)
    p_detail.add_glyph(source, glyph)
    p_detail.add_glyph(source, rects_detail)

    p_detail.grid.visible = False
    p_detail.xaxis.major_label_text_font_style = "bold"
    p_detail.yaxis.minor_tick_line_width = 0
    p_detail.yaxis.major_tick_line_width = 0

    p = gridplot([[p_head],[p_detail]], toolbar_location='below')

    @pn.depends(coloring_mode_btn.param.value, watch=True)
    def update_color(coloring_mode_btn):
        global current_colors_mode

        current_colors_mode = coloring_mode_btn
        colors = get_colors(seqs, coloring_mode_btn)
        source.data.update(colors=colors)

    return p

def get_colors(seqs, mode=1):
    """Make colors for bases in sequence. Color 'C' and 'G' green if they form 'CG'."""
    text = [i for s in list(seqs) for i in s]

    if mode == 1:
        clrs = {'A':'white','T':'white','G':'white','C':'white','-':'white'}
        colors = [clrs[i] for i in text]
        # Re-color 'C' and 'G' in 'CG' to green
        for i in range(len(colors) - 1):
            if text[i] == 'C' and text[i+1] == 'G':
                colors[i] = 'green'   # 'C' in 'CG' is green
                colors[i+1] = 'green'  # 'G' in 'CG' is green
        return colors
    elif mode == 2:
        clrs =  {'A':'red','T':'green','G':'orange','C':'blue','-':'white'}
        colors = [clrs[i] for i in text]
        return colors



w1 = pnw.RadioBoxGroup(name='Gap method', options={'Define minimum gap': 1, 'Set target average gap (optimize min gap)': 2})
w2 = pnw.IntInput(start=0, placeholder='Your value', sizing_mode='stretch_width')

w3_start = pnw.IntInput(start=0, placeholder='Initial length', sizing_mode='stretch_width')
w3_end = pnw.IntInput(start=0, placeholder='Terminal length', sizing_mode='stretch_width')
# Define labels
label_x = pn.pane.HTML('Initial', styles={'width': '100%', 'text-align': 'center'})
label_y = pn.pane.HTML('Terminal', styles={'width': '100%', 'text-align': 'center'})

# Define GridBox layout
w3_packaging_signals = pn.GridBox(
    label_x, label_y,
    w3_start, w3_end,
    ncols=2, sizing_mode='stretch_width'
)

new_seq_id = pnw.TextInput(placeholder='Enter an identifier', sizing_mode='stretch_width')
mutate_btn = pnw.Button(name='Mutate', button_type='primary')
A_rich_toggle = pnw.RadioButtonGroup(name='A-rich', options={'Yes': 1, 'No': 2}, value=2, button_type='default', button_style='solid')

modifiers = pn.FlexBox('# Mutation settings', '## CpG gap options', w1, w2, '## Protect terminal nucleotides', w3_packaging_signals, '## Increase A-richness in remaining codons', A_rich_toggle, '## Name the mutated sequence', new_seq_id, pn.FlexBox(mutate_btn, styles={'margin-top': '10px'}, justify_content='center'), flex_direction='column', styles={'flex-wrap': 'wrap'})
modifiers.visible = False

footer = pn.pane.HTML("""
            Aleksander Sułkowski 2023
            """, styles={'width': '100%', 'text-align': 'center', 'margin-top': '20px', 'margin-bottom': '20px'})

template.sidebar.extend([modifiers, footer])

bokeh_pane = pn.pane.Bokeh(sizing_mode='stretch_width', margin=10)

sequences = []
def download_alignment():

    output = StringIO()
    SeqIO.write(sequences, output, "fasta")
    fasta_str = output.getvalue()
    fasta_bytes = fasta_str.encode()

    return BytesIO(fasta_bytes)


global fileDownloadAdded
fileDownloadAdded = False
download_btn = pnw.FileDownload(callback=download_alignment, filename='alignment.fasta', button_type='primary', label='Download alignment')

def test():
    p1 = figure(height = 400, width = 400, title = 'My Line Plot')
    p1.line([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], line_width = 2, color = 'navy', alpha = 0.5)
    
    p2 = figure(height=200, width=200)
    p2.line([1, 2, 3, 4, 5], [5, 2, 3, 7, 6], line_width=2, color='red', alpha=0.5)

    p3 = figure(height=200, width=200)
    p3.line([1, 2, 3, 4, 5], [2, 5, 8, 2, 7], line_width=2, color='green', alpha=0.5)

    p4 = figure(height=200, width=200)
    p4.line([1, 2, 3, 4, 5], [7, 2, 4, 9, 1], line_width=2, color='black', alpha=0.5)

    # Make a grid
    grid = gridplot([[p1, p2], [p3, p4]])
    return grid

global current_colors_mode
current_colors_mode = 1

successful_load_dummy = pnw.Checkbox(value=False, visible=False)

file_input = pnw.FileInput(accept='.fasta,.aln')
fasta_bytes = file_input.param.value
# Define function to load FASTA file and display sequences
@pn.depends(fasta_bytes)
def load_fasta(fasta_bytes):
    global sequences, current_coloring_mode, fileDownloadAdded, df

    sequences = []
    # Reset the dataframe
    df = pd.DataFrame(columns=['Name', 'Average CpG gap', 'Mutation settings', 'CpG count', 'CpG change', 'A change'])
    df.index.name = 'Index'

    # Get the loaded file bytes
    if fasta_bytes:
        # Decode the bytes to a string
        fasta_string = fasta_bytes.decode('utf-8')
        # Create a StringIO object from the string
        fasta_io = StringIO(fasta_string)
        # Parse the sequences
        sequences = list(SeqIO.parse(fasta_io, "fasta"))
        
        seqs.value = sequences
        coloring_mode_btn.value = current_coloring_mode
        
        sequence_lengths = [len(sequence) for sequence in sequences]
        if len(sequences) > 1 and len(set(sequence_lengths)) > 1:

            top_message.visible = True
            return

        for sequence in sequences:
            sequence.seq = sequence.seq.upper()

         # Add data to the DataFrame for each sequence
        for sequence in sequences:
            try:
                gene = Gene(sequence.seq, gap_method=None)
            except:
                pn.state.notifications.warning('Your file contains invalid codons.')
                top_message.visible = True
                return
            else:
                update_table(average_CpG_gap=gene.original_average_gap, sequence_id=sequence.id, CpG_count=gene.count_CpGs(gene.original_sequence))


        # Display the sequences
        p = view_alignment()
        bokeh_pane.object = p
        
        modifiers.visible = True
        top_message.visible = False
        accessory.append(coloring_mode_btn)

        download_btn.filename = f"{sequences[0].id}-recoding-aln.fasta"
        if len(sequences) > 1 and not fileDownloadAdded:
            accessory.append(download_btn)
            fileDownloadAdded = True
        elif len(sequences) <= 1 and fileDownloadAdded:
            accessory.remove(download_btn)
            fileDownloadAdded = False
        
        successful_load_dummy.value = not successful_load_dummy.value     

       
def mutate(event, gap_method, option_value, protection_start_value, protection_end_value, identifier, A_rich):
    global sequences
    global fileDownloadAdded

    original_id = sequences[0].id
    
    if gap_method == 1:
        gene = Gene(sequences[0].seq, gap_method, minimum_CpG_gap=option_value, packaging_signal_length_beginning=protection_start_value, packaging_signal_length_end=protection_end_value)
        description = f"Recoded {original_id} sequence to increase CpG's with a minimum gap of {option_value}"
        if A_rich == True:
            gene.mutate_A_rich()
            description = f"Recoded {original_id} sequence to increase CpG's with a minimum gap of {option_value}. Mutated remaining codons to A-rich versions, if possible."
        sequences.append(SeqRecord(gene.new_sequence, id=identifier, name=identifier, description=description))
    elif gap_method == 2:
        gene = Gene(sequences[0].seq, gap_method, minimum_CpG_gap=None, desired_CpG_gap=option_value, packaging_signal_length_beginning=protection_start_value, packaging_signal_length_end=protection_end_value)
        resultant_avg_gap = round(gene.calculate_average_gap())
        description = f"Recoded {original_id} sequence to increase CpG's at an average gap of {resultant_avg_gap}"
        if A_rich == True:
            gene.mutate_A_rich()
            description = f"Recoded {original_id} sequence to increase CpG's at an average gap of {resultant_avg_gap}. Mutated remaining codons to A-rich versions, if possible."
        sequences.append(SeqRecord(gene.new_sequence, id=identifier, name=identifier, description=description))


    bokeh_pane.object = view_alignment() #sequences, colors_mode=current_colors_mode

    if len(sequences) > 1 and not fileDownloadAdded:
        accessory.append(download_btn)
        fileDownloadAdded = True
    elif len(sequences) <= 1 and fileDownloadAdded:
        accessory.remove(download_btn)
        fileDownloadAdded = False
    
    if gap_method == 1:
        mutation_settings = f"Minimum CpG gap set as {option_value}. "
    elif gap_method == 2:
        mutation_settings = f"Desired average CpG gap set as {option_value}.\\nAlgorithm-optimized minimal gap value is {gene.closest_gap}."
    if protection_start_value > 0 and protection_end_value == 0:
        mutation_settings = mutation_settings + f"\\nProtected {protection_start_value} initial nucleotides. "
    elif protection_start_value > 0 and protection_end_value > 0:
        mutation_settings = mutation_settings + f"\\nProtected {protection_start_value} initial and {protection_end_value} terminal nucleotides. "
    elif protection_start_value == 0 and protection_end_value > 0:
        mutation_settings = mutation_settings + f"\\nProtected {protection_start_value} terminal nucleotides. "
    if A_rich:
        mutation_settings = mutation_settings + f"\\nA-enriched remaining codons. "
    
    CpG_change = gene.calculate_CpG_abundance_change()
    if CpG_change > 0:
        CpG_change = f"+{CpG_change}"
    A_change = gene.calculate_A_abundance_change()
    if A_change > 0:
        A_change = f"+{A_change}"

    update_table(average_CpG_gap=gene.calculate_average_gap(), sequence_id=identifier, CpG_count=gene.count_CpGs(gene.new_sequence), mutation_settings=mutation_settings.strip(), CpG_abundance_change=f"{CpG_change}%", A_abundance_change=f"{A_change}%")

def on_mutate_btn_click(event):
    gap_method = w1.value
    option_value = w2.value

    if w2.value == None:
        pn.state.notifications.warning('Please enter the mutation value.')
        return
    
    if w3_start.value == None:
        w3_start.value = 0
    if w3_end.value == None:
        w3_end.value = 0
    protection_start_value = w3_start.value
    protection_end_value = w3_end.value
    identifier = new_seq_id.value.strip()

    if A_rich_toggle.value == 1:
        A_rich = True
    elif A_rich_toggle.value == 2:
        A_rich = False

    if identifier is None or identifier == '':
        pn.state.notifications.warning('Please enter at least one character.')
    else:
        found = False
        for seq in sequences:
            if seq.id == identifier:
                found = True
                pn.state.notifications.warning('A sequence with this name already exists.')
                return

        if not found:
            mutate(event=event, gap_method=gap_method, option_value=option_value, protection_start_value=protection_start_value, protection_end_value=protection_end_value, identifier=identifier, A_rich=A_rich)


top = pn.FlexBox('## Load a new FASTA File', file_input, successful_load_dummy, justify_content='center')
top_message = pn.pane.Markdown("""
                           Please load an in-frame sequence, of a length divisible by 3.

                           If loading a Fasta file with more than one sequence, ensure that sequence lengths are equal.

                           In such case, the sequence at the top of the file will be used in the recoding algorithm.
                           """)
accessory = pn.FlexBox(justify_content='space-between')

visualization = pn.FlexBox(top, top_message, load_fasta, bokeh_pane, accessory, flex_direction='column', align_items='center', sizing_mode='stretch_width')

def update_table(average_CpG_gap, sequence_id, CpG_count, mutation_settings=None, CpG_abundance_change=None, A_abundance_change=None):
    global df
    
    data = []

    data = pd.DataFrame([{
        'Name': sequence_id,
        'Average CpG gap': average_CpG_gap,
        'Mutation settings': mutation_settings,
        'CpG count': CpG_count,
        'CpG change': CpG_abundance_change,
        'A change': A_abundance_change,
    }])
    
    # Concatenate new data to existing DataFrame
    df = pd.concat([df, data], ignore_index=True)
    df.index.name = 'Index'
    df.fillna('NA', inplace=True)
    tabulator.value = df
    
df = pd.DataFrame(columns=['Name', 'Average CpG gap', 'Mutation settings', 'CpG count', 'CpG change', 'A change'])
df.index.name = 'Index'
formatter = {
    'Mutation settings': {'type': 'textarea'}
}
tabulator = pnw.Tabulator(df, formatters=formatter, layout='fit_data', disabled=True)

global tableDownloadAdded
tableDownloadAdded = False

table = pn.FlexBox(tabulator, flex_direction='column', align_content='center', sizing_mode='stretch_width')


template.main.extend([
    visualization,
    table,
    ]
)

successful_load_dummy.jscallback(value="""
    if (document.getElementById('sidebar-button') &&
    document.getElementById('sidebar').classList.contains('hidden')) {
        document.getElementById('sidebar-button').click();
    }
    """)

mutate_btn.param.watch(on_mutate_btn_click, 'clicks')

template.servable()

await write_doc()
  `

  try {
    const [docs_json, render_items, root_ids] = await self.pyodide.runPythonAsync(code)
    self.postMessage({
      type: 'render',
      docs_json: docs_json,
      render_items: render_items,
      root_ids: root_ids
    })
  } catch(e) {
    const traceback = `${e}`
    const tblines = traceback.split('\n')
    self.postMessage({
      type: 'status',
      msg: tblines[tblines.length-2]
    });
    throw e
  }
}

self.onmessage = async (event) => {
  const msg = event.data
  if (msg.type === 'rendered') {
    self.pyodide.runPythonAsync(`
    from panel.io.state import state
    from panel.io.pyodide import _link_docs_worker

    _link_docs_worker(state.curdoc, sendPatch, setter='js')
    `)
  } else if (msg.type === 'patch') {
    self.pyodide.globals.set('patch', msg.patch)
    self.pyodide.runPythonAsync(`
    state.curdoc.apply_json_patch(patch.to_py(), setter='js')
    `)
    self.postMessage({type: 'idle'})
  } else if (msg.type === 'location') {
    self.pyodide.globals.set('location', msg.location)
    self.pyodide.runPythonAsync(`
    import json
    from panel.io.state import state
    from panel.util import edit_readonly
    if state.location:
        loc_data = json.loads(location)
        with edit_readonly(state.location):
            state.location.param.update({
                k: v for k, v in loc_data.items() if k in state.location.param
            })
    `)
  }
}

startApplication()