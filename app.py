import objects
import random, io
import string
import numpy as np

from Bio.Seq import Seq
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
from bokeh.io import show

import panel as pn
import panel.widgets as pnw

css = """
:root {
  --accent-foreground-cut-rest: #107bba;
}
"""
pn.extension(raw_css=[css], notifications=True)


def view_alignment(aln, fontsize="9pt", plot_width=800):
    """Bokeh sequence alignment view"""

    #make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]
    text = [i for s in seqs for i in s]
    colors = get_colors(seqs)    
    N = len(seqs[0])
    S = len(seqs)    
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    #use recty for rect coords with an offset
    recty = gy+.5
    h= 1/S
    #now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs)*15+50
    x_range = Range1d(0,N+1, bounds='auto')
    if N>100:
        viewlen=100
    else:
        viewlen=N
    #view_range is for the close up view
    view_range = Range1d(0, viewlen, bounds=(0, N+1))
    tools="xpan, xwheel_zoom, reset"

    #entire sequence view (no text, with zoom)
    p_head = figure(title=None, width= plot_width, height=50,
               x_range=x_range, y_range=(0,S), tools=tools,
               min_border=0, toolbar_location='below', active_scroll='xwheel_zoom')
    rects_head = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                 line_color=None, fill_alpha=0.6)
    p_head.add_glyph(source, rects_head)
    p_head.yaxis.visible = False
    p_head.grid.visible = False  

    #sequence text view with ability to scroll along x axis
    p_detail = figure(title=None, width=plot_width, height=plot_height,
                x_range=view_range, y_range=ids, tools="xpan, reset, xwheel_pan",
                min_border=0, toolbar_location='below', active_scroll='xwheel_pan', lod_threshold=None)          
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                text_font_size=fontsize)
    rects_detail = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)
    p_detail.add_glyph(source, glyph)
    p_detail.add_glyph(source, rects_detail)

    p_detail.grid.visible = False
    p_detail.xaxis.major_label_text_font_style = "bold"
    p_detail.yaxis.minor_tick_line_width = 0
    p_detail.yaxis.major_tick_line_width = 0

    p = gridplot([[p_head],[p_detail]], toolbar_location='below')
    return p

def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i for s in list(seqs) for i in s]
    clrs =  {'A':'red','T':'green','G':'orange','C':'blue','-':'white'}
    colors = [clrs[i] for i in text]
    return colors

# Create a FileInput widget
file_input = pnw.FileInput()
w1 = pnw.RadioBoxGroup(name='Gap method', options={'Define minimum gap': 1, 'Define average gap': 2})
w2 = pnw.IntInput(start=0, placeholder='Your value', sizing_mode='stretch_width')

w3_start = pnw.IntInput(start=0, value=0, placeholder='Initial length', sizing_mode='stretch_width')
w3_end = pnw.IntInput(start=0, value=0, placeholder='Terminal length', sizing_mode='stretch_width')
w3_packaging_signals = pn.Row(w3_start, w3_end)

new_seq_id = pnw.TextInput(placeholder='Enter an identifier')
mutate_btn = pnw.Button(name='Mutate',width=100,button_type='primary')

modifiers = pn.Column('## Mutation options', w1, w2, '## Define packaging signals', w3_packaging_signals, '## Name mutated sequence', new_seq_id, mutate_btn)
modifiers.visible = False

template = pn.template.FastListTemplate(
    title="SynRecoder",
    sidebar=['## Load a FASTA File', file_input, modifiers
],
    theme_toggle=False,
    accent='#A01346',
)


bokeh_pane = pn.pane.Bokeh(height=300,margin=10)

global fileDownloadAdded
fileDownloadAdded = False

sequences = []
# Define function to load FASTA file and display sequences
def load_fasta(event):
    global sequences
    # Get the loaded file path
    fasta_file = file_input.filename
    if fasta_file:
        # Parse the sequences
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        # Display the sequences
        p = view_alignment(sequences)
        bokeh_pane.object = p
        
        modifiers.visible = True
        main.visible = True


def mutate(event, gap_method, option_value, protection_start_value, protection_end_value, identifier):
    global sequences
    global fileDownloadAdded

    original_id = sequences[0].id
    
    if gap_method == 1:
        gene = objects.Gene(sequences[0].seq, gap_method, minimum_CpG_gap=option_value, packaging_signal_length_beginning=protection_start_value, packaging_signal_length_end=protection_end_value)
        sequences.append(SeqRecord(gene.new_sequence, id=identifier, name=identifier, description=f"Recoded {original_id} sequence with a minimum gap of {option_value}"))
    elif gap_method == 2:
        gene = objects.Gene(sequences[0].seq, gap_method, minimum_CpG_gap=None, desired_CpG_gap=option_value, packaging_signal_length_beginning=protection_start_value, packaging_signal_length_end=protection_end_value)
        resultant_avg_gap = round(gene.calculate_average_gap())
        sequences.append(SeqRecord(gene.new_sequence, id=f"{original_id}-avggap-{resultant_avg_gap}", name=f"{original_id}-avggap-{resultant_avg_gap}", description=f"Recoded {original_id} sequence with an average gap of {resultant_avg_gap}"))


    bokeh_pane.object = view_alignment(sequences)

    if not fileDownloadAdded:
            main.append(pn.widgets.FileDownload(callback=download_alignment, filename=f"{original_id}-recoding-aln.fasta", button_style='solid', button_type='primary', label='Download alignment', align='end'))
            fileDownloadAdded = True


def on_mutate_btn_click(event):
    gap_method = w1.value
    option_value = w2.value
    protection_start_value = w3_start.value
    protection_end_value = w3_end.value
    identifier = new_seq_id.value.strip()

    if identifier is None or identifier == '':
        pn.state.notifications.warning('Please enter at least one character.')
    else:
        found = False
        for seq in sequences:
            if seq.id == identifier:
                found = True
                pn.state.notifications.warning('A sequence with this name already exists.')
                break

        if not found:
            mutate(event=event, gap_method=gap_method, option_value=option_value, protection_start_value=protection_start_value, protection_end_value=protection_end_value, identifier=identifier)

def download_alignment():
    filename = f"{sequences[0].id}-recoding-aln.fasta"
    SeqIO.write(sequences, filename, "fasta")
    return filename


main = pn.Column(bokeh_pane, sizing_mode='stretch_height')

template.main.extend([
    main
    ]
)
main.visible = False

# Add event listener to the FileInput widget
file_input.param.watch(load_fasta, 'value')
mutate_btn.param.watch(on_mutate_btn_click, 'clicks')

template.servable()
