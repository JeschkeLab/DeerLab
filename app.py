
# Libraries for basic OS operations
import base64
import os
from urllib.parse import quote as urlquote

# Libraries for the Dash web GUI
from flask import Flask, send_from_directory
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output,State

# Libraries for PlotLy graphics in Dash
import plotly.express as px
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go

# Libraries for DeerLab functionality
import numpy as np
import time
import math
from deerlab import *
from deerload import deerload
import matplotlib.pyplot as plt

# temporary server directory to upload files
UPLOAD_DIRECTORY = "./app_uploaded_files"
if not os.path.exists(UPLOAD_DIRECTORY):
    os.makedirs(UPLOAD_DIRECTORY)

# CSS formatting
external_stylesheets = ['http://dash-gallery.plotly.host/dash-manufacture-spc-dashboard/assets/base-styles.css?m=1583852344.0']

#=====================================================================================
#=========================      BROWSER APP GUI                =======================
#=====================================================================================
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Prepare empty containers for the graphics
fig = make_subplots(rows=1, cols=2)
tab = go.Figure()

app.layout = html.Div(children=[
    dcc.Store(id='experiment_data',storage_type='session',data={'t':[],'Vr':[],'Vi':[]}),
    html.Center(children=[
        html.Img(src=app.get_asset_url('logo_dark.png'),width=400),

        html.Div(children='A simpe of example of a DeerLab-based web dashboard'),
        html.Div(
            [
                html.H4("Upload"),
                dcc.Upload(
                    id="upload-data",
                    children=html.Div(
                        ["Drag and drop or click to select a file to upload."]
                    ),
                    style={
                        "width": "100%",
                        "height": "60px",
                        "lineHeight": "60px",
                        "borderWidth": "1px",
                        "borderStyle": "dashed",
                        "borderRadius": "5px",
                        "textAlign": "center",
                        "margin": "10px",
                    },
                    multiple=True,
                ),
                html.Ul(id="file-list"),
            ],
            style={"max-width": "500px"},
        ),
        html.H4("Model"),
        html.Div("Define you dipolar signal model:"),
        html.Div(
            [
                dcc.Dropdown(id='dd_model_dropdown',
                    options=[
                        {'label': 'Tikhonov', 'value':  'tikh'},
                        {'label': 'dd_gauss', 'value':  'dd_gauss'},
                        {'label': 'dd_gauss2', 'value': 'dd_gauss2'}
                    ],
                    value='tikh',
                    clearable=False,
                    style={"max-width": "400px"},
                ), 
                dcc.Dropdown(id='bg_model_dropdown',
                    options=[
                        {'label': 'bg_hom3d',       'value': 'bg_hom3d'},
                        {'label': 'bg_hoomfractal', 'value': 'bg_hoomfractal'},
                        {'label': 'bg_exp',         'value': 'bg_exp'}
                    ],
                    value='bg_hom3d',
                    clearable=False,
                    style={"max-width": "400px"},
                ), 
                dcc.Dropdown(id='ex_model_dropdown',
                    options=[
                        {'label': '4pDEER',         'value': 'ex_4pdeer'},
                        {'label': '5pDEER',         'value': 'ex_5pdeer'},
                        {'label': '4pDEER (w/ 2+1)',         'value': 'ex_ovl4pdeer'}
                    ],
                    value='ex_4pdeer',
                    clearable=False,
                    style={"max-width": "400px"},
                ), 
            ],
            className='row'
        ),
        html.Button('Submit', id='submit-val', n_clicks=0),
    dcc.Graph(id='example-graph',figure=fig),

    dcc.Graph(id='example-table',figure=tab)
    ]),

])
#=====================================================================================
#=====================================================================================
#=====================================================================================
#=====================================================================================


@app.callback(
    [Output('example-graph', 'figure'),Output('example-table', 'figure')],
    [Input('submit-val', 'n_clicks')],
    [State('experiment_data', 'data'),State('dd_model_dropdown', 'value')])
def update_plots(n_clicks, data, dd_model):
#=====================================================================================
    t = data['t']
    Vr = data['Vr']
    if not Vr:
        fig = make_subplots(rows=1, cols=2)
        tab = go.Figure()
    else:
        fig,tab = processDEER(t,Vr,dd_model)
    return fig,tab
#=====================================================================================

def save_file(name, content):
#=====================================================================================
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(os.path.join(UPLOAD_DIRECTORY, name), "wb") as fp:
        fp.write(base64.decodebytes(data))
#=====================================================================================

def delete_file(name):
#=====================================================================================
    """Delete a file uploaded with Plotly Dash."""
    fullfilename = os.path.join(UPLOAD_DIRECTORY, name) 
    os.remove(fullfilename)
#=====================================================================================


def uploaded_files():
#=====================================================================================
    """List the files in the upload directory."""
    files = []
    for filename in os.listdir(UPLOAD_DIRECTORY):
        path = os.path.join(UPLOAD_DIRECTORY, filename)
        if os.path.isfile(path):
            files.append(filename)
    return files
#=====================================================================================



@app.callback(
    [Output("file-list", "children"), Output("experiment_data","data")],
    [Input("upload-data", "filename"), Input("upload-data", "contents")],
    [State("experiment_data","data")]
)
def upload_external_data(uploaded_filenames, uploaded_file_contents,data):
#=====================================================================================
    """Upload and parse data from external files"""

    # Upload all files to the server and write
    if uploaded_filenames is not None and uploaded_file_contents is not None:
        for name, content in zip(uploaded_filenames, uploaded_file_contents):
            save_file(name, content)

    # Get a list of the files in the server
    files = uploaded_files()

    # Parse files
    for name in files:
        if 'DSC' in name:
            # laod the data from the BRUKER files on the server
            t,V,parameters = deerload(os.path.join(UPLOAD_DIRECTORY, name))
            # package the real/imaginary parts separately to since complex numpy arrays are not JSON-convertible
            data = {'t':t,'Vr':np.real(V),'Vi':np.real(V)}
        # Delete the file from the server
        delete_file(name)

    # If no files are uploaded return just a message
    if len(files) == 0:
        return ([html.Li("No files yet!")],data)
    else:
        items = [html.Li(filename) for filename in files]
        return items,data
#=====================================================================================


def processDEER(t,V,dd_model):
    def Kmodel(p,t,r):

        # Unpack parameters
        lam = p[0]
        k = p[1]

        # Get background
        B = bg_exp(t,k)

        # Generate 4pDEER kernel
        dr = r[2] - r[1]
        K = dipolarkernel(t,r,lam,B)
        return K


    t = np.squeeze(np.asarray(t))
    V = np.squeeze(np.asarray(V))

    t = t/1000

    # phase correction
    V = np.real(V)
    # zerotime correction
    t = t - t[np.argmax(V)]

    r = np.linspace(2,6,100)

    #--------------------------
    # Non-linear parameters:
    #--------------------------
    #       lam  c0
    #--------------------------
    par0 = [0.5, 0.5 ] # Start values
    lb   = [ 0 , 0.05] # lower bounds
    ub   = [ 1 ,  1  ] # upper bounds

    #--------------------------
    # Linear parameters: 
    #--------------------------
    #          Pfit
    #--------------------------
    lbl = np.zeros(len(r)) # Non-negativity constraint of P
    ubl = []

    Amodel = lambda p: Kmodel(p,t,r)

    # Run SNLLS optimization
    parfit, Pfit, uq = snlls(V,Amodel,par0,lb,ub,lbl,ubl)

    # Get fitted model
    Vfit = Kmodel(parfit,t,r)@Pfit

    V = V/np.max(Vfit)
    Vfit = Vfit/np.max(Vfit)
    dr = np.mean(np.diff(r))
    Pnorm = np.sum(Pfit)/dr
    Pfit = Pfit/Pnorm
    P95 = uq.ci(95,'lin')/Pnorm
    P50 = uq.ci(50,'lin')/Pnorm

    fig = make_subplots(rows=1, cols=2)


    # -----------------------------------------------------
    # Time-domain plot
    # -----------------------------------------------------
    fig.add_trace(
        go.Scatter(x=t,y=V,mode="markers",name="Vexp"),
        row=1, col=1
    )

    fig.add_trace(
        go.Scatter(x=t,y=Vfit,mode="lines",name="Vfit"),
        row=1, col=1
    )

    # edit axis labels
    fig['layout']['xaxis']['title']='time [us]'
    fig['layout']['yaxis']['title']='V(t)'
    # -----------------------------------------------------

    # -----------------------------------------------------
    # Distance-domain plot
    # -----------------------------------------------------

    # Fitted distribution
    fig.add_trace(
        go.Scatter(x=r,y=Pfit,mode="lines",name="Pfit",line_color='blue'),
        row=1, col=2
    )

    # Uncertainty 95%-ci
    fig.add_trace(go.Scatter(x=r,y=P95[:,0],mode="lines",name="95%-CI",line_color='blue',line_width=0),
        row=1, col=2
    )
    fig.add_trace(
        go.Scatter(x=r,y=P95[:,1],
                mode="lines",name="95%-CI",fill='tonexty',line_color='blue',line_width=0),
        row=1, col=2
    )

    # Uncertainty 50%-ci
    fig.add_trace(go.Scatter(x=r,y=P50[:,0],mode="lines",name="50%-CI",line_color='blue',line_width=0),
        row=1, col=2
    )
    fig.add_trace(
        go.Scatter(x=r,y=P50[:,1],
                mode="lines",name="50%-CI",fill='tonexty',line_color='blue',line_width=0),
        row=1, col=2
    )

    # edit axis labels
    fig['layout']['xaxis2']['title']='distance [nm]'
    fig['layout']['yaxis2']['title']='P(r)'

    fig.update_layout(showlegend=False)

    parci = uq.ci(95,'nonlin')
    paramnames = ['Modulation depth','Background decay rate']
    tab = go.Figure(data=[go.Table(
        header=dict(values=['Parameter','Fit','95% upper','95% lower'],
                    align='left'),
        cells=dict(values=[paramnames, np.round(parfit,3), np.round(parci[:,0],3), np.round(parci[:,1],3)],
                align='left'))
    ])
    tab.update_layout(width=1000, height=300)


    return fig,tab

if __name__ == "__main__":
    app.run_server(debug=True)


