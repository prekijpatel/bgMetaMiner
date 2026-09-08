from dash import dcc, html, Dash, Input, Output, ctx
import dash_bootstrap_components as dbc
import pandas as pd
import numpy as np
from dash_graphs import Choropleth_map, Assembly_level_bar, Annotation_bar, Submission_year_line, Sequencing_technologies_scatter, Coverage_bar, ANI_scatter, Total_genes_hist, CDSs_hist, Non_coding_hist, Pseudogenes_hist, Isolation_source_treemap, N50L50_scatter
import logging
import dash_daq as daq
import os
import webbrowser
from datetime import datetime
import sys
import logging

# Redirecting standard output to dummy file handles to avoid AttributeError

logging.basicConfig(filename="dash.log", level=logging.DEBUG, format='%(asctime)s:%(levelname)s - %(message)s')


def meta_mined(df: pd.DataFrame, saving_file_path:str = None):
    logging.info(f"Dash started!")

    if saving_file_path is None:
        saving_file_path = os.getcwd()
        logging.debug(f"Since the `saving file path` is `None`, the current working directory is choosen. The WD is {saving_file_path}")
    else:
        logging.debug(f"the current saving path: {saving_file_path}")

    # for country drop-down
    df_for_cdo = df[pd.notnull(df['country_common_name'])]

    country_dropdown_options = df_for_cdo['country_common_name'].dropna().astype(str).unique()

    country_dropdown_options.sort()
    country_dropdown_options = np.insert(country_dropdown_options, 0, 'World')

    # for sequencing technology drop-down
    sequencing_technology_options = df['Categorized_sequencing_technologies'].dropna().astype(str).unique()
    sequencing_technology_options.sort()
    sequencing_technology_options = np.insert(sequencing_technology_options, 0, 'All')

    # rename the column if source_y is present
    if 'source_y' in df.columns:
        logging.debug(f"`source_y` column is present in the dataframe. Renaming it to `source`")
        df.rename(columns={'source_y': 'source'}, inplace=True)
    elif 'source' in df.columns:
        logging.debug(f"`source` column only is there in the dataframe.")

    # for sliders making sure the columns are numeric
    df['Coverage_Depth'] = pd.to_numeric(df['Coverage_Depth'], errors='coerce') # making sure this column is numeric

    df['Total_genes'] = pd.to_numeric(df['Total_genes'], errors='coerce')

    df['Protein-coding_genes'] = pd.to_numeric(df['Protein-coding_genes'], errors='coerce')

    df['Non-coding_genes'] = pd.to_numeric(df['Non-coding_genes'], errors='coerce')

    df['Protein-coding_genes'] = pd.to_numeric(df['Protein-coding_genes'], errors='coerce')

    df['Pseudogenes'] = pd.to_numeric(df['Pseudogenes'], errors='coerce')
    
    df['ANI_best_match_score'] = pd.to_numeric(df['ANI_best_match_score'], errors='coerce')

    df['ANI_best_matched_assembly\'s_coverage'] = pd.to_numeric(df['ANI_best_matched_assembly\'s_coverage'], errors='coerce')

    df['Contig_N50'] = pd.to_numeric(df['Contig_N50'], errors='coerce')

    df['Contig_L50'] = pd.to_numeric(df['Contig_L50'], errors='coerce')

    flask_logger = logging.getLogger('werkzeug')
    flask_logger.setLevel(logging.ERROR)
    
    # print("Initializing app...")
    meta_mined = Dash(__name__, external_stylesheets=[dbc.themes.LUX])

    meta_mined.layout = dcc.Loading(
        id = 'loading-fullpage',
        type='circle',
        fullscreen=True,
        overlay_style={"visibility": "visible"},
        delay_hide=800,
        delay_show=800,
        children=[
            html.Div(
                children=[
                    html.H1(
                        'Meta-mined!',
                        style={'textAlign': 'center', 'fontSize': '2.5em', 'color': '#333', 'margin': '20px 0'}
                    ),
                    html.Hr(style={'border': '1px solid #ccc', 'margin': '10px 0'}),
                    
                    # LED Display for genome count
                    html.Div(
                        style={'display': 'flex', 'justifyContent': 'center', 'marginBottom': '20px'},
                        children=[
                            daq.LEDDisplay(
                                id='genome-count',
                                value=len(df),
                                color="#5e5eff",
                                backgroundColor="#f9f9f9",
                                size=30,
                            ),
                        ],
                    ),
                    
                    # First row: Radio buttons for atypical and suppressed assembly options
                    html.Div(
                        style={'display': 'flex', 'gap': '20px', 'marginBottom': '20px'},
                        children=[
                            html.Div(
                                style={
                                    'flex': 1, 'padding': '10px', 'backgroundColor': '#f9f9f9', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    dcc.RadioItems(
                                        id='atypical-radio',
                                        options=[
                                            {'label': 'Include Atypical Assemblies', 'value': 'all'},
                                            {'label': 'Exclude Atypical Assemblies', 'value': 'no_atypical'},
                                            {'label': 'Only Atypical Assemblies', 'value': 'only_atypical'},
                                        ],
                                        value='all',
                                        labelStyle={'display': 'block', 'margin': '5px 0'},
                                    ),
                                    html.Div(
                                        id='atypical-radio-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555'},
                                    ),
                                ],
                            ),
                            html.Div(
                                style={
                                    'flex': 1, 'padding': '10px', 'backgroundColor': '#f9f9f9', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    dcc.RadioItems(
                                        id='suppressed-radio',
                                        options=[
                                            {'label': 'Include Suppressed Assemblies', 'value': 'all'},
                                            {'label': 'Exclude Suppressed Assemblies', 'value': 'no_suppressed'},
                                            {'label': 'Only Suppressed Assemblies', 'value': 'only_suppressed'},
                                        ],
                                        value='all',
                                        labelStyle={'display': 'block', 'margin': '5px 0'},
                                    ),
                                    html.Div(
                                        id='suppressed-radio-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555'},
                                    ),
                                ],
                            ),
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '10px', 'backgroundColor': '#f9f9f9', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Assembly Submission Year', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='submission-year-slider',
                                        min=1980,
                                        max=2025,
                                        step=1,
                                        marks={i: str(i) for i in range(1980, 2026, 5)},
                                        value=[1980, 2025],
                                    ),
                                    html.Div(
                                        id='submission-year-slider-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                ],
                            ),
                        ],
                    ),
                    
                    # Second row: Choropleth map for isolate's geographical location
                    html.Div(
                        style={'display': 'flex', 'gap': '20px', 'marginBottom': '20px'},
                        children=[
                            html.Div(
                                style={
                                    'flex': 3, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Geographical Distribution of Isolates', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Dropdown(
                                        id='country-dropdown',
                                        options=country_dropdown_options,
                                        value='World',
                                        placeholder='Select a country',
                                        style={'marginBottom': '5px', 'justifyContent': 'center'}
                                    ),
                                    html.Div(
                                        style={'display': 'flex', 'justifyContent': 'center'},
                                        children=[
                                            dcc.Graph(
                                                id='choropleth-map',
                                                figure=Choropleth_map(df=df, selected_country='World'),
                                            )
                                        ]
                                    ),
                                ],
                            ),
                        ],
                    ),

                    # Third row: Graph for assembly level and annotations
                    html.Div(
                        style={'display': 'flex', 'gap': '20px', 'marginBottom': '20px'},
                        children=[
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Assembly Level:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Checklist(
                                        id='assembly-level-checklist',
                                        options=[
                                            {'label': 'Complete Genome', 'value': 'Complete Genome'},
                                            {'label': 'Chromosome', 'value': 'Chromosome'},
                                            {'label': 'Scaffold', 'value': 'Scaffold'},
                                            {'label': 'Contig', 'value': 'Contig'},
                                        ],
                                        value=['Complete Genome', 'Chromosome', 'Scaffold', 'Contig'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                    ),
                                    html.Div(
                                        style={'display': 'flex', 'justifyContent': 'center'},
                                        children=[
                                            dcc.Graph(
                                                id='assembly-level-graph',
                                                figure=Assembly_level_bar(
                                                    df=df,
                                                    selected_assembly_levels=['Complete Genome', 'Chromosome', 'Scaffold', 'Contig']
                                                ),
                                            ),
                                        ]
                                    ),
                                ],
                            ),
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Annotation From:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Checklist(
                                        id='annotation-checklist',
                                        options=[
                                            {'label': 'GenBank', 'value': 'GenBank'},
                                            {'label': 'NCBI RefSeq', 'value': 'NCBI RefSeq'},
                                            {'label': 'Others', 'value': 'Others'},
                                            {'label': 'No Annotation', 'value': 'No Annotation'},
                                        ],
                                        value=['GenBank', 'NCBI RefSeq', 'Others', 'No Annotation'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                    ),
                                    html.Div(
                                        style={'display': 'flex', 'justifyContent': 'center'},
                                        children=[
                                            dcc.Graph(
                                                id='annotation-graph',
                                                figure=Annotation_bar(
                                                    df=df,
                                                    show_annotations_from=['GenBank', 'NCBI RefSeq', 'Others', 'No Annotation']
                                                ),
                                            ),
                                        ],
                                    ),
                                ],
                            ),

                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Submissions over the years:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Graph(
                                        id='submission-year-scatter',
                                        figure=Submission_year_line(df=df),
                                    ),
                                ],
                            ),
                        ],
                    ),
                    
                    # Fourth row: Sequencing technologies and coverage depth
                    html.Div(
                        style={'display': 'flex', 'gap': '20px', 'marginBottom': '20px'},
                        children=[
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Sequencing Technologies', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Dropdown(
                                        id='sequencing-technology-dropdown',
                                        options=[{'label': i, 'value': i} for i in sequencing_technology_options],
                                        multi=True,
                                        value=['All'],
                                        placeholder='Select a sequencing technology',
                                        style={'marginBottom': '5px'}
                                    ),
                                    dcc.Graph(
                                        id='sequencing-technologies-line',
                                        figure=Sequencing_technologies_scatter(df=df, selected_sequencing_technologies=['All']),
                                    ),
                                ],
                            ),
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Coverage Depth', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='coverage-slider',
                                        min=0,
                                        max=5000,
                                        step=25,
                                        marks={i: str(i) for i in range(0, 5001, 500)},
                                        value=[0, 5000],
                                    ),
                                    dcc.Checklist(
                                        id='coverage-checklist',
                                        options=[
                                            {'label': 'Include coverage depth of >5000', 'value': '> 5000'},
                                        ],
                                        value=['> 5000'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='coverage-slider-output1',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Checklist(
                                        id='coverage-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no coverage data', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='coverage-slider-output2',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    html.Div(
                                        style={'display': 'flex', 'justifyContent': 'center'},
                                        children=[
                                            dcc.Graph(
                                                id='coverage-depth-bar',
                                                figure=Coverage_bar(df=df),
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                        ],
                    ),

                    # Fifth row: Average Nucleotide Identity and Contig L50 - N50
                    html.Div(
                        style={'display': 'flex', 'gap': '20px', 'marginBottom': '20px'},
                        children=[
                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Average Nucleotide Identity', style={'marginBottom': '10px', 'color': '#333'}),
                                    html.H6('% Identity', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='ani-identity-slider',
                                        min=0,
                                        max=100,
                                        step=1,
                                        marks={i: str(i) for i in range(0, 101, 10)},
                                        value=[0, 100],
                                    ),
                                    dcc.Checklist(
                                        id='ani-identity-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no ANI % Identity data', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='ani-identity-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    html.H6('% Coverage', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='ani-coverage-slider',
                                        min=0,
                                        max=100,
                                        step=1,
                                        marks={i: str(i) for i in range(0, 101, 10)},
                                        value=[0, 100],
                                    ),
                                    dcc.Checklist(
                                        id='ani-coverage-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no ANI %coverage data', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='ani-coverage-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    html.Div(
                                        id='ani-sliders-outputs',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Graph(
                                        id='ani-scatter',
                                        figure=ANI_scatter(df=df),
                                    ),
                                ],
                            ),

                            html.Div(
                                style={
                                    'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                                    'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                                },
                                children=[
                                    html.H5('Contig L50 - N50:', style={'marginBottom': '10px', 'color': '#333'}),
                                    html.H6('Contig N50', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='n50-slider',
                                        min=int(df['Contig_N50'].min()),
                                        max=int(df['Contig_N50'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Contig_N50'].max()), 1000000)},
                                        value=[0, int(df['Contig_N50'].max())],
                                    ),
                                    html.Br(),
                                    html.H6('Contig L50', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.RangeSlider(
                                        id='l50-slider',
                                        min=int(df['Contig_L50'].min()),
                                        max=int(df['Contig_L50'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Contig_L50'].max()), 500)},
                                        value=[0, int(df['Contig_L50'].max())],
                                    ),
                                    html.Div(
                                        id='n50l50-slider-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    html.Br(),
                                    dcc.Graph(
                                        id='n50l50-scatter',
                                        figure=N50L50_scatter(df=df),
                                    ),
                                ],
                            ),
                        ],
                    ),

                    # Sixth row: histograms for total genes, CDSs, non-coding genes, and pseudogenes
                    html.Div(
                        style={
                            'display': 'grid',
                            'gridTemplateColumns': 'repeat(4, 1fr)',  
                            'gap': '20px',  
                            'marginBottom': '20px',
                            # # 'flex': 2, 'padding': '15px', 'backgroundColor': '#ffffff', 
                            # # 'borderRadius': '8px', 'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)'
                            'padding': '20px',
                            'backgroundColor': '#f8f9fa',
                            'borderRadius': '10px'
                        },
                        children=[
                            
                            html.H5('Annotations:', style={'gridColumn': 'span 4', 'marginBottom': '10px', 'color': '#333', 'marginLeft':'20px'}),
                            # html.Br(),

                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    dcc.RangeSlider(
                                        id = 'total-gene-slider',
                                        min = int(df['Total_genes'].min()),
                                        max = int(df['Total_genes'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Total_genes'].max()), 2000)},
                                        value=[0, int(df['Total_genes'].max())],
                                    ),
                                    html.Div(
                                        id='total-gene-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Checklist(
                                        id='total-gene-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no total gene count', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='total-gene-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Graph(
                                        id='total-gene-hist',
                                        figure=Total_genes_hist(df=df),
                                    ),
                                ]
                            ),

                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    dcc.RangeSlider(
                                        id='cds-slider',
                                        min=int(df['Protein-coding_genes'].min()),
                                        max=int(df['Protein-coding_genes'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Protein-coding_genes'].max()), 2000)},
                                        value=[0, int(df['Protein-coding_genes'].max())],
                                    ),
                                    html.Div(
                                        id='cds-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Checklist(
                                        id='cds-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no CDSs count', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='cds-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Graph(
                                        id='cds-hist',
                                        figure=CDSs_hist(df=df),
                                    ),
                                ]
                            ),

                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    dcc.RangeSlider(
                                        id='non-coding-slider',
                                        min=int(df['Non-coding_genes'].min()),
                                        max=int(df['Non-coding_genes'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Non-coding_genes'].max()), 50)},
                                        value=[0, int(df['Non-coding_genes'].max())],
                                    ),
                                    html.Div(
                                        id='non-coding-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Checklist(
                                        id='non-coding-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no non-coding gene count', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='non-coding-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Graph(
                                        id='non-coding-hist',
                                        figure=Non_coding_hist(df=df),
                                    ),
                                ]
                            ),

                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    dcc.RangeSlider(
                                        id='pseudogene-slider',
                                        min=int(df['Pseudogenes'].min()),
                                        max=int(df['Pseudogenes'].max()),
                                        step=1,
                                        marks={i: str(i) for i in range(0, int(df['Pseudogenes'].max()), 500)},
                                        value=[0, int(df['Pseudogenes'].max())],
                                    ),
                                    html.Div(
                                        id='pseudogene-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Checklist(
                                        id='pseudogene-include-null-checklist',
                                        options=[
                                            {'label': 'Include genomes with no Pseudogene count', 'value': 'all'},
                                        ],
                                        value=['all'],
                                        labelStyle={'display': 'inline-block', 'marginRight': '5px', 'fontSize': '0.95em'},
                                        style={'textAlign':'center'}
                                    ),
                                    html.Div(
                                        id='pseudogene-null-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                    dcc.Graph(
                                        id='pseudogene-hist',
                                        figure=Pseudogenes_hist(df=df),
                                    ),
                                ]
                            ),
                        ]
                    ),


                    html.Div(
                        style={
                            'display': 'block',
                            'gap': '20px',
                            'marginBottom': '20px',
                            'padding': '20px',
                        },
                        children=[
                            html.H5('Filter BioProject:', style={'marginBottom': '10px', 'color': '#333',}),
                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    html.H6('Search BioProject by keywords:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Input(
                                        id='bioproject-input',
                                        type='text',
                                        placeholder='Search Keywords separated by commas...',
                                        debounce=True,
                                        style={
                                            'width': '700px',
                                            'height': '40px',         
                                            'fontSize': '16px',       
                                            'padding': '5px 10px',    
                                            'borderRadius': '5px',
                                            'marginBottom': '15px',     
                                        }
                                    ),
                                    html.H6('Selected BioProjects:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Dropdown(
                                        id='bioproject-dropdown',
                                        options=[{'label': d, 'value': d} for d in df['Bioproject_title'].unique()],
                                        placeholder='Filtered BioProjects will appear here...',
                                        searchable=True,
                                        multi=True,
                                    ),
                                    html.Div(
                                        id='bioproject-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                ],
                            ),
                        ],
                    ),

                    html.Div(
                        style={
                            'display': 'block',
                            'gap': '20px',
                            'marginBottom': '20px',
                            'padding': '20px',
                        },
                        children=[
                            html.H5('Filter BioSample:', style={'marginBottom': '10px', 'color': '#333', }),
                            html.Div(
                                style={
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '8px',
                                    'boxShadow': '0 2px 4px rgba(0, 0, 0, 0.1)', 
                                },
                                children=[
                                    html.H6('Search BioSample by keywords:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Input(
                                        id='biosample-input',
                                        type='text',
                                        placeholder='Search Keywords separated by commas...',
                                        debounce=True,
                                        style={
                                            'width': '700px',
                                            'height': '40px',         
                                            'fontSize': '16px',       
                                            'padding': '5px 10px',    
                                            'borderRadius': '5px',
                                            'marginBottom': '15px'    
                                        }
                                    ),
                                    html.H6('Selected BioSamples:', style={'marginBottom': '10px', 'color': '#333'}),
                                    dcc.Dropdown(
                                        id='biosample-dropdown',
                                        options=[{'label': d, 'value': d} for d in df['Biosample_title'].unique()],
                                        placeholder='Filtered BioSamples will appear here...',
                                        searchable=True,
                                        multi=True,
                                    ),
                                    html.Div(
                                        id='biosample-output',
                                        style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'},
                                    ),
                                ],
                            ),
                        ],
                    ),

                    # Eighth row: Isolation source treemap and dropdowns
                    html.Div(
                        style={
                            'display': 'flex',
                            'gap': '20px',
                            'marginBottom': '20px',
                            'padding': '20px',
                            'backgroundColor': '#f8f9fa',
                            'borderRadius': '10px'
                        },
                        children=[
                            # Left Panel - Dropdowns in Two Sections
                            html.Div(
                                style={
                                    'flex': 2,
                                    'padding': '20px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '10px',
                                    'boxShadow': '0 4px 8px rgba(0, 0, 0, 0.1)',
                                    'display': 'flex',
                                    'flexDirection': 'column',
                                    'gap': '15px'
                                },
                                children=[
                                    html.H5('Isolation Source', style={'marginBottom': '10px', 'color': '#333', 'fontWeight': 'bold'}),

                                    # All Dropdowns in One Column First
                                    html.Div(
                                        style={'display': 'flex', 'flexDirection': 'column', 'gap': '15px'},
                                        children=[
                                            html.Div(
                                                children=[
                                                    html.H6('Select Host', style={'marginBottom': '8px', 'color': '#444'}),
                                                    dcc.Dropdown(
                                                        id='identified-host-dropdown',
                                                        options=[
                                                            {'label': 'Hospital-associated', 'value': 'Hospital-associated'},
                                                            {'label': 'Animal-associated', 'value': 'Animal-associated'},
                                                            {'label': 'Environment-associated', 'value': 'Environment-associated'},
                                                            {'label': 'Laboratory-based', 'value': 'Laboratory-based'},
                                                            {'label': 'Unknown', 'value': 'Unknown'}
                                                        ],
                                                        multi=True,
                                                        value=['Hospital-associated', 'Animal-associated', 'Environment-associated', 'Laboratory-based','Unknown'],
                                                        placeholder="Select Host",
                                                    )
                                                ]
                                            ),
                                            html.Div(
                                                children=[
                                                    html.H6('Select Source/Disease Category', style={'marginBottom': '8px', 'color': '#444'}),
                                                    dcc.Dropdown(
                                                        id='source-category-dropdown',
                                                        multi=True,
                                                        placeholder="Select Source/Disease Category",
                                                    )
                                                ]
                                            ),
                                            html.Div(
                                                children=[
                                                    html.H6('Select Source/Disease', style={'marginBottom': '8px', 'color': '#444'}),
                                                    dcc.Dropdown(
                                                        id='source-dropdown',
                                                        multi=True,
                                                        placeholder="Select Source/Disease",
                                                    )
                                                ]
                                            ),
                                            html.Div(
                                                children=[
                                                    html.H6('Select Sample', style={'marginBottom': '8px', 'color': '#444'}),
                                                    dcc.Dropdown(
                                                        id='sample-dropdown',
                                                        multi=True,
                                                        placeholder="Select Sample",
                                                    )
                                                ]
                                            ),
                                        ]
                                    ),
                                ],
                            ),
                            # Right Panel - Treemap Output on Top, Graph Below
                            html.Div(
                                style={
                                    'flex': 1,
                                    'display': 'flex',
                                    'flexDirection': 'column',
                                    'alignItems': 'center',
                                    'padding': '15px',
                                    'backgroundColor': '#ffffff',
                                    'borderRadius': '10px',
                                    'boxShadow': '0 4px 8px rgba(0, 0, 0, 0.1)',
                                    'gap': '15px'
                                },
                                children=[
                                    # Treemap Output (Text or Additional Info)
                                    html.Div(
                                        id='treemap-output',
                                        style={
                                            'textAlign': 'center',
                                            'fontSize': '0.9em',
                                            'color': '#555',
                                            'padding': '10px',
                                            'backgroundColor': '#f1f1f1',
                                            'borderRadius': '8px',
                                            'width': '100%'
                                        },
                                    ),

                                    # Treemap Graph
                                    dcc.Graph(
                                        id='isolation-source-treemap',
                                        figure=Isolation_source_treemap(df=df),
                                        style={'width': '100%', 'height': '100%'}
                                    ),
                                ]
                            ),
                        ]
                    ),

        
                    html.Div(
                        children=[
                            html.Div(
                                children=[
                                    html.Button(
                                        "Save data",
                                        id='save-button',
                                        n_clicks=0,
                                        style={
                                            'backgroundColor': '#007bff',  
                                            'color': 'white', 
                                            'border': 'none', 
                                            'borderRadius': '5px', 
                                            'padding': '10px 20px', 
                                            'fontSize': '1em', 
                                            'cursor': 'pointer',
                                            'boxShadow': '0px 4px 6px rgba(0, 0, 0, 0.1)'
                                        }
                                    ),
                                ],
                                style={
                                    'display': 'flex',  
                                    'justifyContent': 'center',  
                                    'alignItems': 'center',  
                                    'height': '100px'  
                                }
                            ),
                            html.Div(
                                id='save-button-output',
                                style={'textAlign': 'center', 'fontSize': '0.9em', 'color': '#555', 'marginTop': '10px'}
                            ),
                        ]
                    ),
                ],
            ),
        ],
    )
    # Define callback
    @meta_mined.callback(
        [Output('choropleth-map', 'figure'),
        Output('assembly-level-graph', 'figure'),
        Output('annotation-graph', 'figure'),
        Output('submission-year-slider-output', 'children'),
        Output('atypical-radio-output', 'children'),
        Output('suppressed-radio-output', 'children'),
        Output('genome-count', 'value'),
        Output('sequencing-technologies-line', 'figure'),
        Output('submission-year-scatter', 'figure'),
        Output('coverage-depth-bar', 'figure'),
        Output('coverage-slider-output1', 'children'),
        Output('coverage-slider-output2', 'children'),
        Output('ani-scatter', 'figure'),
        Output('ani-identity-null-output', 'children'),
        Output('ani-sliders-outputs', 'children'),
        Output('ani-coverage-null-output', 'children'),
        Output('n50l50-slider-output', 'children'),
        Output('n50l50-scatter', 'figure'), 
        [Output('total-gene-output', 'children'),
         Output('total-gene-hist', 'figure')],
        Output('total-gene-null-output', 'children'),
        [Output('cds-output', 'children'),
         Output('cds-hist', 'figure')],
        Output('cds-null-output', 'children'),
        [Output('non-coding-output', 'children'),
         Output('non-coding-hist', 'figure')],
        Output('non-coding-null-output', 'children'),
        [Output('pseudogene-output', 'children'),
         Output('pseudogene-hist', 'figure')],
        Output('pseudogene-null-output', 'children'),
        [Output('bioproject-dropdown', 'options'),
        Output('bioproject-dropdown', 'value')],
        Output('bioproject-output', 'children'),
        [Output('biosample-dropdown', 'options'),
        Output('biosample-dropdown', 'value')],
        Output('biosample-output', 'children'),
        Output('source-category-dropdown', 'options'),
        Output('source-dropdown', 'options'),
        Output('sample-dropdown', 'options'),
        [Output('treemap-output', 'children'),
        Output('isolation-source-treemap', 'figure')],
        Output('save-button-output', 'children'),],
        [Input('country-dropdown', 'value'),
        Input('assembly-level-checklist', 'value'),
        Input('annotation-checklist', 'value'),
        Input('submission-year-slider', 'value'),
        Input('atypical-radio', 'value'),
        Input('suppressed-radio', 'value'),
        Input('genome-count', 'id'),
        Input('sequencing-technology-dropdown', 'value'),
        Input('submission-year-scatter', 'id'),
        Input('coverage-slider', 'value'),
        Input('coverage-checklist', 'value'),
        Input('coverage-include-null-checklist', 'value'),
        Input('ani-identity-slider', 'value'),
        Input('ani-identity-include-null-checklist', 'value'),
        Input('ani-coverage-slider', 'value'),
        Input('ani-coverage-include-null-checklist', 'value'),
        Input('n50-slider', 'value'),
        Input('l50-slider', 'value'),
        Input('total-gene-slider', 'value'),
        Input('total-gene-include-null-checklist', 'value'),
        Input('cds-slider', 'value'),
        Input('cds-include-null-checklist', 'value'),
        Input('non-coding-slider', 'value'),
        Input('non-coding-include-null-checklist', 'value'),
        Input('pseudogene-slider', 'value'),
        Input('pseudogene-include-null-checklist', 'value'),
        Input('bioproject-input', 'value'),
        Input('bioproject-dropdown', 'value'),
        Input('biosample-input', 'value'),
        Input('biosample-dropdown', 'value'),
        Input('identified-host-dropdown', 'value'),
        Input('source-category-dropdown', 'value'),
        Input('source-dropdown', 'value'),
        Input('sample-dropdown', 'value'),
        Input('save-button', 'n_clicks')]
    )

    # print("we are at update_dash")
    def update_dash(selected_country, selected_assembly_level, show_annotations_from, year_range, atypical_radio, suppressed_radio, genome_count_id, selected_sequencing_technologies, submission_year_line_id, coverage_range, coverage_checklist, coverage_include_null_checklist, ani_identity_range, ani_identity_include_null_checklist, ani_coverage_range, ani_coverage_include_null_checklist, n50_range, l50_range,total_genes_range, total_genes_include_null_checklist, cds_range, cds_include_null_checklist, non_coding_range, non_coding_include_null_checklist, pseudogene_range, pseudogene_include_null_checklist, bioproject_input, bioproject_dropdown_values, biosample_input, biosample_dropdown_values, identified_host, source_category, source, sample, n_save_button_clicks):
        modified_df = df.copy()

        logging.debug(f"Shape of the modified_df before filtering: {modified_df.shape}")

        for col in modified_df.columns:
            if modified_df[col].apply(lambda x: isinstance(x, dict)).any():
                logging.debug(f"Converting {col} data to string values as dictionaries may cause problem in filtering/removing duplicates.")
                modified_df[col] = modified_df[col].apply(str)


        if selected_country == 'World':
            modified_df = modified_df
            logging.debug(f"Seleted country is World. Shape of the modified_df after filtering by country: {modified_df.shape}")
        else:
            modified_df = modified_df[modified_df['country_common_name'] == selected_country]
            logging.debug(f"Selected country is {selected_country}. Shape of the modified_df after filtering by country: {modified_df.shape}")

        modified_df = modified_df[modified_df['Assembly_level'].isin(selected_assembly_level)]
        logging.debug(f"Selected assembly levels: {selected_assembly_level}")
        logging.debug(f"Shape of the modified_df after filtering by assembly level: {modified_df.shape}")

        modified_df = modified_df[modified_df['Annotation_category'].isin(show_annotations_from)]
        logging.debug(f"Selected annotation categories: {show_annotations_from}")
        logging.debug(f"Shape of the modified_df after filtering by annotation category: {modified_df.shape}")

        
        modified_df = modified_df[(modified_df['Submission_year'] >= year_range[0]) & (modified_df['Submission_year'] <= year_range[1])]


        if atypical_radio == 'no_atypical':
            earlier_length = modified_df.shape[0]
            modified_df = modified_df[modified_df['Assmbly_atypical?'] == 'No']
            atypical_output = f"Number of atypical assemblies excluded: {earlier_length - modified_df.shape[0]}"

        elif atypical_radio == 'only_atypical':
            earlier_length = modified_df.shape[0]
            modified_df = modified_df[modified_df['Assmbly_atypical?'] == 'Yes']
            atypical_output = f"Number of atypical assemblies: {modified_df.shape[0]}"
        else:
            modified_df = modified_df
            atypical_output = None

        if suppressed_radio == 'no_suppressed':
            all_length = modified_df.shape[0]
            modified_df = modified_df[modified_df['Assembly_status'] == 'current']
            suppressed_output = f"Number of suppressed assemblies excluded: {all_length - modified_df.shape[0]}"
        elif suppressed_radio == 'only_suppressed':
            all_length = modified_df.shape[0]
            modified_df = modified_df[modified_df['Assembly_status'] == 'suppressed']
            suppressed_output = f"Number of suppressed assemblies: {modified_df.shape[0]}"
        else:
            modified_df = modified_df
            suppressed_output = None

        if 'All' in selected_sequencing_technologies:
            modified_df = modified_df
            logging.debug(f"Selected sequencing technology: {selected_sequencing_technologies}")
            logging.debug(f"Shape of the modified_df after filtering by sequencing technology: {modified_df.shape}")
        else:
            modified_df = modified_df[modified_df['Categorized_sequencing_technologies'].isin(selected_sequencing_technologies)]
            logging.debug(f"Selected sequencing technology: {selected_sequencing_technologies}")
            logging.debug(f"Shape of the modified_df after filtering by sequencing technology: {modified_df.shape}")
        
        logging.debug(f"Shape of the modified_df after sorting sequencing technologies: {modified_df.shape}")
        coverage_null_text=""
        if coverage_range[0] == 0 & coverage_range[1] == 5000:            
            if '> 5000' in coverage_checklist and 'all' in coverage_include_null_checklist:
                modified_df = modified_df  
            elif '> 5000' in coverage_checklist and 'all' not in coverage_include_null_checklist:
                coverage_null_text = f"Genomes with no coverage data are excluded."
                coverage_df_above_5000 = modified_df[modified_df['Coverage_Depth'] > 5000]
                coverage_df_below_5000 = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
                modified_df = pd.concat([coverage_df_above_5000, coverage_df_below_5000], ignore_index=True)
            elif '> 5000' not in coverage_checklist and 'all' in coverage_include_null_checklist:
                coverage_null_df = modified_df[pd.isnull(modified_df['Coverage_Depth'])]
                coverage_df_below_5000 = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
                modified_df = pd.concat([coverage_null_df, coverage_df_below_5000], ignore_index=True)
            elif '> 5000' not in coverage_checklist and 'all' not in coverage_include_null_checklist:
                coverage_null_text = f"Genomes with no coverage data are excluded."
                modified_df = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
        else:
            if '> 5000' in coverage_checklist and 'all' in coverage_include_null_checklist:
                coverage_null_df = modified_df[pd.isnull(modified_df['Coverage_Depth'])]
                coverage_df_above_5000 = modified_df[modified_df['Coverage_Depth'] > 5000]
                coverage_df_below_5000 = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
                modified_df = pd.concat([coverage_null_df, coverage_df_above_5000, coverage_df_below_5000], ignore_index=True)
            elif '> 5000' in coverage_checklist and 'all' not in coverage_include_null_checklist:
                coverage_null_text = f"Genomes with no coverage data are excluded."
                coverage_df_above_5000 = modified_df[modified_df['Coverage_Depth'] > 5000]
                coverage_df_below_5000 = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
                modified_df = pd.concat([coverage_df_above_5000, coverage_df_below_5000], ignore_index=True)
            elif '> 5000' not in coverage_checklist and 'all' in coverage_include_null_checklist:
                coverage_null_df = modified_df[pd.isnull(modified_df['Coverage_Depth'])]
                coverage_df_below_5000 = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]
                modified_df = pd.concat([coverage_null_df, coverage_df_below_5000], ignore_index=True)
            elif '> 5000' not in coverage_checklist and 'all' not in coverage_include_null_checklist:
                coverage_null_text = f"Genomes with no coverage data are excluded."
                modified_df = modified_df[(modified_df['Coverage_Depth'] >= coverage_range[0]) & (modified_df['Coverage_Depth'] <= coverage_range[1])]

        for col in modified_df.columns:
            if modified_df[col].apply(lambda x: isinstance(x, dict)).any():
                logging.debug(f"Column '{col}' contains dictionaries.")

        logging.debug(f"Shape of the modified_df after coverage filters: {modified_df.shape}")
        ani_identity_null_text=""
        if ani_identity_range[0] == 0 & ani_identity_range[1] == 100:
            if 'all' in ani_identity_include_null_checklist:
                modified_df = modified_df
            else:
                ani_identity_null_text = f"{len(modified_df[pd.isnull(modified_df['ANI_best_match_score'])])} Genomes with no % ANI identity data are excluded."
                modified_df = modified_df[pd.notnull(modified_df['ANI_best_match_score'])]
        else:
            if 'all' in ani_identity_include_null_checklist:
                ANI_identity_null_df = modified_df[pd.isnull(modified_df['ANI_best_match_score'])]
                selected_ANI_identity_df = modified_df[(modified_df['ANI_best_match_score'] >= ani_identity_range[0]) & (modified_df['ANI_best_match_score'] <= ani_identity_range[1])]
                modified_df = pd.concat([selected_ANI_identity_df, ANI_identity_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                ani_identity_null_text = f"{len(modified_df[pd.isnull(modified_df['ANI_best_match_score'])])} Genomes with no % ANI identity data are excluded."
                modified_df = modified_df[(modified_df['ANI_best_match_score'] >= ani_identity_range[0]) & (modified_df['ANI_best_match_score'] <= ani_identity_range[1])]
        
        ani_coverage_null_text = ""
        if ani_coverage_range[0] == 0 & ani_coverage_range[1] == 100:
            if 'all' in ani_coverage_include_null_checklist:
                modified_df = modified_df                
            else:
                ani_coverage_null_count = len(modified_df[pd.isnull(modified_df['ANI_best_matched_assembly\'s_coverage'])])
                ani_coverage_null_text = f"{ani_coverage_null_count} genomes with no % coverage data are excluded."
                modified_df = modified_df[pd.notnull(modified_df['ANI_best_matched_assembly\'s_coverage'])]
        else:
            if 'all' in ani_coverage_include_null_checklist:
                ANI_coverage_null_df = modified_df[pd.isnull(modified_df['ANI_best_matched_assembly\'s_coverage'])]
                selected_ANI_coverage_df = modified_df[(modified_df['ANI_best_matched_assembly\'s_coverage'] >= ani_coverage_range[0]) & (modified_df['ANI_best_matched_assembly\'s_coverage'] <= ani_coverage_range[1])]
                modified_df = pd.concat([selected_ANI_coverage_df, ANI_coverage_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                ani_coverage_null_count = len(modified_df[pd.isnull(modified_df['ANI_best_matched_assembly\'s_coverage'])])
                ani_coverage_null_text = f"{ani_coverage_null_count} genomes with no % coverage data are excluded."
                modified_df = modified_df[(modified_df['ANI_best_matched_assembly\'s_coverage'] >= ani_coverage_range[0]) & (modified_df['ANI_best_matched_assembly\'s_coverage'] <= ani_coverage_range[1])]

        if n50_range[0] == 0 & n50_range[1] == int(df['Contig_N50'].max()):
            modified_df = modified_df
        else:
            modified_df = modified_df[(modified_df['Contig_N50'] >= n50_range[0]) & (modified_df['Contig_N50'] <= n50_range[1])]
        
        if l50_range[0] == 0 & l50_range[1] == int(df['Contig_L50'].max()):
            modified_df = modified_df
        else:
            modified_df = modified_df[(modified_df['Contig_L50'] >= l50_range[0]) & (modified_df['Contig_L50'] <= l50_range[1])]

        logging.debug(f"Shape of the modified_df after ani filters: {modified_df.shape}")
        total_gene_null_text = ""
        if total_genes_range[0] == 0 & total_genes_range[1] == int(df['Total_genes'].max()):
            if 'all' in total_genes_include_null_checklist:
                modified_df = modified_df
            else:
                total_gene_null_text = f"{len(modified_df[pd.isnull(modified_df['Total_genes'])])} genomes excluded that didn't have total gene count."
                modified_df = modified_df[pd.notnull(modified_df['Total_genes'])]
        else:
            if 'all' in total_genes_include_null_checklist:
                total_gene_null_df = modified_df[pd.isnull(modified_df['Total_genes'])]
                selected_total_gene_df = modified_df[(modified_df['Total_genes'] >= total_genes_range[0]) & (modified_df['Total_genes'] <= total_genes_range[1])]
                modified_df = pd.concat([selected_total_gene_df, total_gene_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                total_gene_null_text = f"{len(modified_df[pd.isnull(modified_df['Total_genes'])])} genomes excluded that didn't have total gene count."
                modified_df = modified_df[(modified_df['Total_genes'] >= total_genes_range[0]) & (modified_df['Total_genes'] <= total_genes_range[1])]

        cds_null_text = ""
        if cds_range[0] == 0 & cds_range[1] == int(df['Protein-coding_genes'].max()):
            if 'all' in cds_include_null_checklist:
                modified_df = modified_df
            else:
                cds_null_text = f"{len(modified_df[pd.isnull(modified_df['Protein-coding_genes'])])} genomes excluded that didn't have CDSs count."
                modified_df = modified_df[pd.notnull(modified_df['Protein-coding_genes'])]
        else:
            if 'all' in cds_include_null_checklist:
                cds_null_df = modified_df[pd.isnull(modified_df['Protein-coding_genes'])]
                selected_cds_df = modified_df[(modified_df['Protein-coding_genes'] >= cds_range[0]) & (modified_df['Protein-coding_genes'] <= cds_range[1])]
                modified_df = pd.concat([selected_cds_df, cds_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                cds_null_text = f"{len(modified_df[pd.isnull(modified_df['Protein-coding_genes'])])} genomes excluded that didn't have CDSs count."
                modified_df = modified_df[(modified_df['Protein-coding_genes'] >= cds_range[0]) & (modified_df['Protein-coding_genes'] <= cds_range[1])]

        non_coding_null_text = ""
        if non_coding_range[0] == 0 & non_coding_range[1] == int(df['Non-coding_genes'].max()):
            if 'all' in non_coding_include_null_checklist:
                modified_df = modified_df
            else:
                non_coding_null_text = f"{len(modified_df[pd.isnull(modified_df['Non-coding_genes'])])} genomes excluded that didn't have non-coding gene count."
                modified_df = modified_df[pd.notnull(modified_df['Non-coding_genes'])]
        else:
            if 'all' in non_coding_include_null_checklist:
                non_coding_null_df = modified_df[pd.isnull(modified_df['Non-coding_genes'])]
                selected_non_coding_df = modified_df[(modified_df['Non-coding_genes'] >= non_coding_range[0]) & (modified_df['Non-coding_genes'] <= non_coding_range[1])]
                modified_df = pd.concat([selected_non_coding_df, non_coding_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                non_coding_null_text = f"{len(modified_df[pd.isnull(modified_df['Non-coding_genes'])])} genomes excluded that didn't have non-coding gene count."
                modified_df = modified_df[(modified_df['Non-coding_genes'] >= non_coding_range[0]) & (modified_df['Non-coding_genes'] <= non_coding_range[1])]
        
        pseudogene_null_text = ""
        if pseudogene_range[0] == 0 & pseudogene_range[1] == int(df['Pseudogenes'].max()):
            if 'all' in pseudogene_include_null_checklist:
                modified_df = modified_df
            else:
                pseudogene_null_text = f"{len(modified_df[pd.isnull(modified_df['Pseudogenes'])])} genomes excluded that didn't have pseudogene count"
                modified_df = modified_df[pd.notnull(modified_df['Pseudogenes'])]
        else:
            if 'all' in pseudogene_include_null_checklist:
                pseudogene_null_df = modified_df[pd.isnull(modified_df['Pseudogenes'])]
                selected_pseudogene_df = modified_df[(modified_df['Pseudogenes'] >= pseudogene_range[0]) & (modified_df['Pseudogenes'] <= pseudogene_range[1])]
                modified_df = pd.concat([selected_pseudogene_df, pseudogene_null_df], ignore_index=True)
                modified_df.drop_duplicates()
            else:
                pseudogene_null_text = f"{len(modified_df[pd.isnull(modified_df['Pseudogenes'])])} genomes excluded that didn't have pseudogene count"
                modified_df = modified_df[(modified_df['Pseudogenes'] >= pseudogene_range[0]) & (modified_df['Pseudogenes'] <= pseudogene_range[1])]
        
        logging.debug(f"Shape of the modified_df after gene count filters: {modified_df.shape}")

        # filtering the dataframe based on keywords in bioproject and biosample names  
        
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0] if ctx.triggered else None

        logging.debug(f"Triggered ID: {triggered_id}")
        
        bioproject_options = [{'label': i, 'value': i} for i in modified_df['Bioproject_title'].unique()]

        if bioproject_input:

            bioproject_list = modified_df['Bioproject_title'].tolist()

            keywords = [word.strip().lower() for word in bioproject_input.split(',') if word.strip()]
            
            logging.debug(f"Input for bioproject: {bioproject_input}")

            list_based_on_bioproject_input = [title for title in bioproject_list if any(kw in title.lower() for kw in keywords)]

        else:
            bioproject_dropdown = []
            list_based_on_bioproject_input = []
            logging.debug(f"No filter was used for Bioproject")

        if triggered_id == 'bioproject-input':
            bioproject_dropdown = list(set(bioproject_dropdown_values or []) | set(list_based_on_bioproject_input))
            logging.debug("User typed in input box — adding matches to current dropdown selection")
        elif triggered_id == 'bioproject-dropdown':
            bioproject_dropdown = list(set(bioproject_dropdown_values)) or []
            logging.debug("User changed dropdown — using dropdown values only")
        else:
            if bioproject_dropdown_values:
                # bioproject_dropdown = list(set(bioproject_dropdown_values)) or []
                bioproject_dropdown = [title for title in modified_df['Bioproject_title'].unique() if title in bioproject_dropdown_values]
                logging.debug("Fallback: using whatever is available in dropdown")

        if bioproject_dropdown != []:
            modified_df = modified_df[modified_df['Bioproject_title'].isin(bioproject_dropdown)]
            logging.debug(f"Filtered {len(modified_df)} rows using: {triggered_id}")
        else:
            modified_df = modified_df
            logging.debug("No filtering applied for bioproject")
        

        biosample_options = [{'label': i, 'value': i} for i in modified_df['Biosample_title'].unique()]

        if biosample_input:
            biosample_list = modified_df['Biosample_title'].tolist()

            keywords_2 = [word.strip().lower() for word in biosample_input.split(',') if word.strip()]
            
            logging.debug(f"Input for biosample: {biosample_input}")

            list_based_on_biosample_input = [title for title in biosample_list if any(kw in title.lower() for kw in keywords_2)]
        
        else:
            biosample_dropdown = []
            list_based_on_biosample_input = []
            logging.debug(f"No filter was used for Biosample")
        
        if triggered_id == 'biosample-input':
            biosample_dropdown = list(set(biosample_dropdown_values or []) | set(list_based_on_biosample_input))
            logging.debug("User typed in input box — adding matches to current dropdown selection")
        elif triggered_id == 'biosample-dropdown':
            biosample_dropdown = list(set(biosample_dropdown_values)) or []
            logging.debug("User changed dropdown — using dropdown values only")
        else:
            if biosample_dropdown_values:
                # biosample_dropdown = list(set(biosample_dropdown_values)) or []
                biosample_dropdown = [title for title in modified_df['Biosample_title'].unique() if title in biosample_dropdown_values]
                logging.debug("Fallback: using whatever is available in dropdown")
        
        if biosample_dropdown != []:
            modified_df = modified_df[modified_df['Biosample_title'].isin(biosample_dropdown)]
            logging.debug(f"Filtered {len(modified_df)} rows using: {triggered_id}")

        else:
            modified_df = modified_df
            logging.debug("No filtering applied for biosample")

        modified_df = modified_df[modified_df['identified_host'].isin(identified_host)]
        logging.debug(f"Selected hosts: {identified_host}")
        logging.debug(f"Shape of the modified_df after filtering by host: {modified_df.shape}")

        source_category_options = [{'label': i, 'value': i} for i in modified_df['source_category'].unique()]

        if source_category:
            modified_df = modified_df[modified_df['source_category'].isin(source_category)]
            logging.debug(f"Selected source categories: {source_category}")
            logging.debug(f"Shape of the modified_df after filtering by source category: {modified_df.shape}")
            source_options = [{'label': i, 'value': i} for i in modified_df['source'].unique()]
        else:
            modified_df = modified_df
            source_options = []
        
        if source:
            modified_df = modified_df[modified_df['source'].isin(source)]
            logging.debug(f"Selected sources: {source}")
            logging.debug(f"Shape of the modified_df after filtering by source: {modified_df.shape}")

            sample_options = [{'label': i, 'value': i} for i in modified_df['sample'].unique()]
        else:
            modified_df = modified_df
            sample_options = []

        if sample:
            modified_df = modified_df[modified_df['sample'].isin(sample)]
            logging.debug(f"Selected samples: {sample}")
            logging.debug(f"Shape of the modified_df after filtering by sample: {modified_df.shape}")
        else:
            modified_df = modified_df

        genome_count = modified_df.shape[0]

        logging.debug(f"The genome count is {genome_count}")

        location_map = Choropleth_map(df=modified_df, selected_country=selected_country)
        
        assembly_level_graph = Assembly_level_bar(df=modified_df, selected_assembly_levels=selected_assembly_level)
        
        annotation_graph = Annotation_bar(df=modified_df, show_annotations_from=show_annotations_from)

        submission_year_text = f"Assembly Submission Year Range: {year_range[0]} - {year_range[1]}"

        sequencing_technologies_over_the_years = Sequencing_technologies_scatter(df=modified_df, selected_sequencing_technologies=selected_sequencing_technologies)

        submissions_over_the_years = Submission_year_line(df=modified_df)

        coverage_depth_bar = Coverage_bar(df=modified_df)

        coverage_depth_text = f"Coverage Depth range: {modified_df['Coverage_Depth'].min()}X to {modified_df['Coverage_Depth'].max()}X"
        
        ani_score_scatter = ANI_scatter(df=modified_df)
        
        ani_score_text = f"Selected % Identity: {ani_identity_range[0]}  to {ani_identity_range[1]}, Selected % Coverage: {ani_coverage_range[0]} to {ani_coverage_range[1]}" 

        n50l50_output_text= f"Contig N50: {n50_range[0]} - {n50_range[1]}, Contig L50: {l50_range[0]} - {l50_range[1]}"

        n50l50_scatter = N50L50_scatter(df=modified_df)
        
        total_gene_dist = Total_genes_hist(df=modified_df)
        total_gene_text = f"Selected 'Total gene' range: {total_genes_range[0]} to {total_genes_range[1]}"

        cds_dist = CDSs_hist(df=modified_df)
        cds_text = f"Selected 'CDSs' range: {cds_range[0]} to {cds_range[1]}"
        
        non_coding_dist = Non_coding_hist(df=modified_df)
        non_coding_text = f"Selected 'Non-coding gene' range: {non_coding_range[0]} to {non_coding_range[1]}"

        pseudogene_dist = Pseudogenes_hist(df=modified_df)
        pseudogene_text = f"Selected 'Pseudogene' range: {pseudogene_range[0]} to {pseudogene_range[1]}"

        bioproject_output = f"{len(modified_df['Bioproject_title'].unique())} BioProjects are selected. They contain a total of {len(modified_df)}assemblies."
        
        biosample_output = f"Selected BioSamples contain a total of {len(modified_df)} assemblies."
        
        treemap_df = modified_df.copy()

        treemap_df = treemap_df[treemap_df['identified_host'] != 'Unknown']

        treemap_text = f"Isolation source is Unknown for {len(modified_df[modified_df['identified_host'] == 'Unknown'])} genomes among all selected genomes."
        isolation_treemap = Isolation_source_treemap(df=treemap_df)

        if n_save_button_clicks > 0:
            logging.debug("Button Pressed!")

            try:
                date_time = datetime.now()
                date_time = str(date_time).replace(" ", "_")
                date_time = date_time.replace(":", "-")
                # print(date_time)
                saving_file_name = "filtered_data_" + date_time + ".tsv"
                # print(saving_file_name)
                
                saving_file_as = str(os.path.join(saving_file_path, saving_file_name)) 
                # print(saving_file_as)
                modified_df.to_csv(saving_file_as, index=False, sep='\t')
                logging.info(f"File saved successfully at {os.path.abspath(saving_file_as)}")
                save_message = f"File saved successfully at {os.path.abspath(saving_file_as)}"
            except Exception as e:
                logging.info(f"Error saving file: {str(e)}")
                save_message =  f"Error saving file: {str(e)}"
            
        else:
            save_message = ""



        return location_map, assembly_level_graph, annotation_graph, submission_year_text, atypical_output, suppressed_output, genome_count, sequencing_technologies_over_the_years, submissions_over_the_years, coverage_depth_bar, coverage_depth_text, coverage_null_text, ani_score_scatter, ani_identity_null_text, ani_score_text, ani_coverage_null_text, n50l50_output_text, n50l50_scatter, [total_gene_text, total_gene_dist], total_gene_null_text, [cds_text, cds_dist], cds_null_text, [non_coding_text, non_coding_dist], non_coding_null_text, [pseudogene_text, pseudogene_dist], pseudogene_null_text, [bioproject_options, bioproject_dropdown], bioproject_output, [biosample_options, biosample_dropdown], biosample_output, source_category_options, source_options, sample_options, [treemap_text, isolation_treemap], save_message

    
        # webbrowser.open("http://localhost:8050")
    
    
    # meta_mined.run_server(debug=True, port=8050)
    meta_mined.run(debug=True, port=8050)  # Turn off reloader if inside Jupyter

# if __name__ == '__main__':
# import pickle

# with open("tmp/df.pkl", "rb") as f:
#     df = pickle.load(f)

# logging.info(f"Data loaded!")

# try:
#     with open("tmp/save_path.pkl", "rb") as f:
#         saving_file_path = pickle.load(f)
# except:
#     logging.debug("The saving file path is none.")
#     saving_file_path = None

# logging.info(f"Saving path loaded.")
# # print("Data loaded!")

# logging.info("Running the dashboard for visulization and exploration.")
# # print("Running the dashboard for visulization and exploration.")
# try:
#     # webbrowser.open("http://localhost:8050")
#     # print("putting metaminer in with df!")
#     logging.debug(f"Putting data into metaminer. The save directory is {saving_file_path}.")
#     meta_mined(df, saving_file_path=saving_file_path)
    
# except Exception as e:
#     logging.error(f"Error running dashboard: {e}", exc_info=True)
#     # print(f"Error running dashboard: {e}")
#     sys.exit(1)


if __name__ == '__main__':
    df = pd.read_csv('./tmp/Allmetadata_2025-04-28_18-31-04.951296.tsv', sep='\t', low_memory=False)

    # print(df.shape)Allmetadata_2025-01-28_19-20-08.193277

    # from  utils2 import *

    # df_2 = all_normalization_operations(df)

    print(type(df))

    meta_mined(df, saving_file_path=None)