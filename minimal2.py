from bokeh.plotting import figure
import panel as pn

pn.extension(notifications=True)

# Create a button
button = pn.widgets.Button(name='Click me', button_type='primary')

def create_figure():
    plot = figure(height = 400, width = 400, title = 'My Line Plot')
    plot.line([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], line_width = 2, color = 'navy', alpha = 0.5)
    return plot

# Callback function to be executed when the button is clicked
def button_click(event):
    # Create a simple line plot
    plot = create_figure()

    # Create a Panel object to display the plot
    panel_obj.object = plot

# Attach the callback function to button's click event
button.on_click(button_click)

panel_obj = pn.pane.Bokeh()
# Create a Panel to show the button
panel = pn.Row(button, panel_obj)

# Serve the Panel
panel.servable()