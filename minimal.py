import panel as pn

pn.extension(notifications=True)

# Create a button
button = pn.widgets.Button(name='Click me', button_type='primary')

# Callback function to be executed when the button is clicked
def button_click(event):
    pn.state.notifications.warning("Warning!")

# Attach the callback function to button's click event
button.on_click(button_click)

# Create a Panel to show the button
panel = pn.Row(button)

# Serve the Panel
panel.servable()