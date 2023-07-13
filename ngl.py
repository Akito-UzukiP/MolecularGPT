import urllib.request as urlreq
from dash import Dash, dcc, html, Input, Output, callback
import dash_bio as dashbio
from dash_bio.utils import xyz_reader

app = Dash(__name__)


data = urlreq.urlopen(
    'https://git.io/speck_methane.xyz'
).read().decode('utf-8')

data = xyz_reader.read_xyz(datapath_or_datastring=data, is_datafile=False)

app.layout = html.Div([
    dcc.Dropdown(
        id='default-speck-preset-views',
        options=[
            {'label': 'Default', 'value': 'default'},
            {'label': 'Ball and stick', 'value': 'stickball'}
        ],
        value='default'
    ),
    dashbio.Speck(
        id='default-speck',
        data=data
    ),
])

@callback(
    Output('default-speck', 'presetView'),
    Input('default-speck-preset-views', 'value')
)
def update_preset_view(preset_name):
    return preset_name


if __name__ == '__main__':
    app.run(debug=True)
