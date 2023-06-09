import objects
import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

import panel as pn
import panel.widgets as pnw


pn.extension('tabulator', design='fast', notifications=True)


def view_alignment(aln, fontsize="9pt", colors_mode=1):
    """Bokeh sequence alignment view"""

    #make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]
    text = [i for s in seqs for i in s]
    colors = get_colors(seqs, colors_mode)    
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
    h= 1/S
    #now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
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
               min_border=0, toolbar_location='below', active_scroll='xwheel_zoom')
    rects_head = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                 line_color=None, fill_alpha=0.6)
    p_head.add_glyph(source, rects_head)
    p_head.yaxis.visible = False
    p_head.grid.visible = False  

    #sequence text view with ability to scroll along x axis
    p_detail = figure(title=None, sizing_mode='scale_width', height=plot_height,
                x_range=view_range, y_range=ids[::-1], tools="xpan, reset, xwheel_pan",
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

    @pn.depends(coloring_mode_btn.param.value, watch=True)
    def update_color(coloring_mode_btn):
        global current_colors_mode

        current_colors_mode = coloring_mode_btn
        colors = get_colors(seqs, coloring_mode_btn)
        source.data.update(colors=colors)

    return p

def get_colors(seqs, mode):
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

template = pn.template.FastListTemplate(
    title="SynRecoder",
    sidebar=[modifiers,
             footer,
],
    theme_toggle=False,
    accent='#A01346',
    collapsed_sidebar=True,
    main_max_width='80vw',
)


bokeh_pane = pn.pane.Bokeh(sizing_mode='stretch_width', margin=10)

sequences = []
def download_alignment(event):
    filename = f"{sequences[0].id}-recoding-aln.fasta"
    SeqIO.write(sequences, filename, "fasta")
    return filename

global fileDownloadAdded
fileDownloadAdded = False
download_btn = pnw.FileDownload(callback=download_alignment, button_style='solid', button_type='primary', label='Download alignment')


global current_colors_mode
current_colors_mode = 1

successful_load_dummy = pnw.Checkbox(value=False, visible=False)

coloring_mode_btn = pnw.RadioButtonGroup(name='Coloring mode', options={"Highlight CpG's": 1, 'Color all nucleotides': 2}, button_style='solid', button_type='light')
# Define function to load FASTA file and display sequences
def load_fasta(event):
    global sequences, current_coloring_mode, fileDownloadAdded, df
    
    sequences = []
    # Reset the dataframe
    df = pd.DataFrame(columns=['Name', 'Average CpG gap', 'Mutation settings', 'CpG count', 'CpG change', 'A change'])
    df.index.name = 'Index'

    # Get the loaded file path
    fasta_file = file_input.filename
    if fasta_file:
        # Parse the sequences
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
        
        
        sequence_lengths = [len(sequence) for sequence in sequences]
        if len(sequences) > 1 and len(set(sequence_lengths)) > 1:
            pn.state.notifications.warning('Please enter sequences of equal length.')
            top_message.visible = True
            return
        
        for sequence in sequences:
            sequence.seq = sequence.seq.upper()

         # Add data to the DataFrame for each sequence
        for sequence in sequences:
            try:
                gene = objects.Gene(sequence.seq, gap_method=None)
            except:
                pn.state.notifications.warning('Your file contains invalid codons.')
                top_message.visible = True
                return
            else:
                update_table(average_CpG_gap=gene.original_average_gap, sequence_id=sequence.id, CpG_count=gene.count_CpGs(gene.original_sequence))


        # Display the sequences
        p = view_alignment(sequences, colors_mode=current_colors_mode)
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
        gene = objects.Gene(sequences[0].seq, gap_method, minimum_CpG_gap=option_value, packaging_signal_length_beginning=protection_start_value, packaging_signal_length_end=protection_end_value)
        description = f"Recoded {original_id} sequence to increase CpG's with a minimum gap of {option_value}"
        if A_rich == True:
            gene.mutate_A_rich()
            description = f"Recoded {original_id} sequence to increase CpG's with a minimum gap of {option_value}. Mutated remaining codons to A-rich versions, if possible."
        sequences.append(SeqRecord(gene.new_sequence, id=identifier, name=identifier, description=description))
    elif gap_method == 2:
        gene = objects.Gene(sequences[0].seq, gap_method, minimum_CpG_gap=None, desired_CpG_gap=option_value, packaging_signal_length_beginning=protection_start_value, packaging_signal_length_end=protection_end_value)
        resultant_avg_gap = round(gene.calculate_average_gap())
        description = f"Recoded {original_id} sequence to increase CpG's at an average gap of {resultant_avg_gap}"
        if A_rich == True:
            gene.mutate_A_rich()
            description = f"Recoded {original_id} sequence to increase CpG's at an average gap of {resultant_avg_gap}. Mutated remaining codons to A-rich versions, if possible."
        sequences.append(SeqRecord(gene.new_sequence, id=identifier, name=identifier, description=description))


    bokeh_pane.object = view_alignment(sequences, colors_mode=current_colors_mode)

    if len(sequences) > 1 and not fileDownloadAdded:
        accessory.append(download_btn)
        fileDownloadAdded = True
    elif len(sequences) <= 1 and fileDownloadAdded:
        accessory.remove(download_btn)
        fileDownloadAdded = False
    
    if gap_method == 1:
        mutation_settings = f"Minimum CpG gap set as {option_value}. "
    elif gap_method == 2:
        mutation_settings = f"Desired average CpG gap set as {option_value}.\nAlgorithm-optimized minimal gap value is {gene.closest_gap}."
    if protection_start_value > 0 and protection_end_value == 0:
        mutation_settings = mutation_settings + f"\nProtected {protection_start_value} initial nucleotides. "
    elif protection_start_value > 0 and protection_end_value > 0:
        mutation_settings = mutation_settings + f"\nProtected {protection_start_value} initial and {protection_end_value} terminal nucleotides. "
    elif protection_start_value == 0 and protection_end_value > 0:
        mutation_settings = mutation_settings + f"\nProtected {protection_start_value} terminal nucleotides. "
    if A_rich:
        mutation_settings = mutation_settings + f"\nA-enriched remaining codons. "
    
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



# Create a FileInput widget
file_input = pnw.FileInput(accept='.fasta,.aln')
top = pn.FlexBox('## Load a new FASTA File', file_input, successful_load_dummy, justify_content='center')
top_message = pn.pane.Markdown("""
                           Please load an in-frame sequence, of a length divisible by 3.

                           If loading a Fasta file with more than one sequence, ensure that sequence lengths are equal.

                           In such case, the sequence at the top of the file will be used in the recoding algorithm.
                           """)
accessory = pn.FlexBox(justify_content='space-between')

visualization = pn.FlexBox(top, top_message, bokeh_pane, accessory, flex_direction='column', align_items='center', sizing_mode='stretch_width')

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


# Add event listener to the FileInput widget
file_input.param.watch(load_fasta, 'value', queued=True)

successful_load_dummy.jscallback(value="""
    if (document.getElementById('sidebar-button') &&
    document.getElementById('sidebar').classList.contains('hidden')) {
        document.getElementById('sidebar-button').click();
    }
    """)

mutate_btn.param.watch(on_mutate_btn_click, 'clicks')

template.servable()
