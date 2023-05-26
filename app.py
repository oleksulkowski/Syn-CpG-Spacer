import objects
import random, io
import string
import numpy as np

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


import panel as pn
import panel.widgets as pnw
pn.extension()

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot


def view_alignment(aln, fontsize="20pt", plot_width=800):
    """Bokeh sequence alignment view"""

    #make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]    
    text = [i for s in list(seqs) for i in s]
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
    x_range = Range1d(start=0, end=N+1, bounds=(0, N+1))
    if N>100:
        viewlen=100
    else:
        viewlen=N
    #view_range is for the close up view
    view_range = Range1d(start=0, end=viewlen, bounds=(0, N+1))

    #entire sequence view (no text, with zoom)
    p = figure(title=None, width= plot_width, height=50,
               x_range=x_range, y_range=(0,S), tools="xpan",
               min_border=0, toolbar_location=None)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                 line_color=None, fill_alpha=0.6)
    p.add_glyph(source, rects)
    p.yaxis.visible = False
    p.grid.visible = False  

    #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, width=plot_width, height=plot_height,
                x_range=view_range, y_range=ids, tools="xpan,reset,xwheel_pan",
                min_border=0, toolbar_location=None, active_scroll="xwheel_pan")#, lod_factor=1)          
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                text_font="monospace",text_font_size=fontsize)
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0

    p = gridplot([[p],[p1]], toolbar_location=None)
    return p

def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i for s in list(seqs) for i in s]
    clrs =  {'A':'white','T':'white','G':'white','C':'blue','-':'white'}
    colors = [clrs[i.upper()] for i in text]
    return colors
    
# Create a Markdown pane for the title
title = pn.pane.Markdown('## Load a FASTA File')

# Create a FileInput widget
file_input = pnw.FileInput()
mutate_btn = pnw.Button(name='Mutate',width=100,button_type='primary')

bokeh_pane = pn.pane.Bokeh(height=300,margin=10)

gene = None
sequences = []
# Define function to load FASTA file and display sequences
def load_fasta(event):
    global gene, sequences
    # Get the loaded file path
    fasta_file = file_input.filename
    if fasta_file:
        # Parse the sequences
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        # Display the sequences
        bokeh_pane.object = view_alignment(sequences, fontsize="7pt", plot_width=600)
        gene = objects.Gene(sequences[0].seq, 1)


def mutate(event):
    global gene, sequences
    sequences.append(SeqRecord(Seq(gene.new_sequence), id="Rec", name="Rec", description="Rec"))
    bokeh_pane.object = view_alignment(sequences, fontsize="7pt", plot_width=600)

# Add event listener to the FileInput widget
file_input.param.watch(load_fasta, 'value')

mutate_btn.param.watch(mutate, 'clicks')


top = pn.Row(file_input, mutate_btn)
bottom = pn.Row(bokeh_pane, sizing_mode='stretch_height')
app = pn.Column(title,top,bottom)
app.servable()
