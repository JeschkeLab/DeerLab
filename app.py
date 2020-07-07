import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from dash.dependencies import Input, Output, State

import numpy as np
import time
import math
from deerlab import *
import matplotlib.pyplot as plt

def simulateDEER(noiselvl):
    def Kmodel(p,t,r):

        # Unpack parameters
        lam = p[0]
        k = p[1]

        # Get background
        B = bg_exp(t,k)

        # Generate 4pDEER kernel
        dr = r[2] - r[1]
        K = dipolarkernel(t,r)/dr
        K = (1-lam + lam*K)*B[:,np.newaxis]*dr
        return K

    np.random.seed(1)

    t = np.linspace(-0.5,5,250)
    r = np.linspace(2,6,100)

    # Generate ground truth and input signal
    P = dd_gauss2(r,[3.5, 0.4, 0.4, 4.5, 0.7, 0.6])
    lam = 0.36
    k = 0.085 #uM
    K = Kmodel([lam,k],t,r)
    V = K@P + whitegaussnoise(t,noiselvl)

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
    P95 = uq.ci(95,'lin')
    fig.add_trace(go.Scatter(x=r,y=P95[:,0],mode="lines",name="95%-CI",line_color='blue',line_width=0),
        row=1, col=2
    )
    fig.add_trace(
        go.Scatter(x=r,y=P95[:,1],
                mode="lines",name="95%-CI",fill='tonexty',line_color='blue',line_width=0),
        row=1, col=2
    )

    # Uncertainty 50%-ci
    P50 = uq.ci(50,'lin')
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



external_stylesheets = ['http://dash-gallery.plotly.host/dash-manufacture-spc-dashboard/assets/base-styles.css?m=1583852344.0']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

fig,tab = simulateDEER(0.01)

app.layout = html.Div(children=[
    html.Center(children=[
        html.Img(src=app.get_asset_url('logo_dark.png'),width=400),

        html.Div(children='A simpe of example of a DeerLab-based web dashboard'),
        html.Div("Specify the noise level to simulate the DEER trace:"),
        html.Div(dcc.Input(id='input-on-submit', value = 0.01, type='number')),
        html.Button('Submit', id='submit-val', n_clicks=0),
    dcc.Graph(id='example-graph',figure=fig),

    dcc.Graph(id='example-table',figure=tab)
    ]),

])


@app.callback(
    [Output('example-graph', 'figure'),Output('example-table', 'figure')],
    [Input('submit-val', 'n_clicks')],
    [State('input-on-submit', 'value')])
def update_output(n_clicks, value):
    fig,tab = simulateDEER(value)
    return fig,tab

if __name__ == '__main__':
    app.run_server(debug=True)

