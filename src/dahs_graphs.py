import plotly.express as px
import pandas as pd
import json
import numpy as np


def remove_country_prefix(country):
    return country.split('-')[1]

def Choropleth_map(df:pd.DataFrame, selected_country:str):
    if selected_country == 'World':    
        country_counts = df.groupby(['country_common_name', 'country_three_lettered_name']).size().reset_index(name='counts')
        country_counts['log_count'] = np.log1p(country_counts['counts'])
        no_counts = len(df) - country_counts['counts'].sum()

        fig = px.choropleth(
            country_counts,
            locations='country_three_lettered_name',
            color='log_count',
            hover_name='country_common_name',
            hover_data='counts',  
            title='Global Counts',
            color_continuous_scale=px.colors.sequential.Aggrnyl
        )

        fig.update_layout(
            title=dict(
                text=f"Total isolates: {len(df)}<br>Isolates without geographical data: {no_counts}",
                xref='paper',
                yref='container',
                y = 0.05,
                x = 0.5,  
                font=dict(size=14, color='#202A44'),

            ),
            showlegend=False,
            
            geo=dict(
                showframe=True,
                showcoastlines=False,
                projection_type='equirectangular'
            ),

            coloraxis_colorbar_title='Log Counts',

            height=800,
            width=1200,
        )

    elif selected_country == 'United States':

        filtered_data = df[df['country_common_name'] == 'United States']

        state_counts = filtered_data.groupby(['state_name', 'state_code']).size().reset_index(name='counts')

        no_count = len(filtered_data) - state_counts['counts'].sum()

        state_counts['counts'] = pd.to_numeric(state_counts['counts'], errors='coerce')

        state_counts['state_code'] = state_counts['state_code'].apply(lambda x: remove_country_prefix(x) if '-' in x else x)

        fig = px.choropleth(
            state_counts,
            locations='state_code',
            scope='usa',
            locationmode='USA-states',
            hover_name='state_name',
            hover_data='counts',  
            color='counts',
            title=f'Counts for {selected_country}',
            color_continuous_scale=px.colors.sequential.Sunset
        )

        fig.update_layout(
            title=dict(
                text=f"Total isolates from {selected_country}: {len(filtered_data)}<br>Isolates without states' info.: {no_count}",
                xref='paper',
                yref='container',
                y = 0.05,
                x = 0.5,  
            ),
            showlegend=False,
            coloraxis_colorbar_title='Counts',

            height=800,
            width=1200,
        )
    
    else:
        
        filtered_data = df[df['country_common_name'] == selected_country]

        # filtered_data.to_csv(f'filtered_data_{selected_country}.csv', index=False)

        state_counts = filtered_data.groupby(['state_name', 'state_code']).size().reset_index(name='counts')

        no_count = len(filtered_data) - state_counts['counts'].sum()

        # state_counts.to_csv(f'{selected_country}_state_counts_before_adding_zeroes.csv', index=False)

        json_name = filtered_data['country_three_lettered_name'].iloc[0] + '_states'

        with open(f'./data/geojsons/{json_name}.json', encoding='utf-8') as f:
            country_geojson = json.load(f)

        geo_states = {state['properties']['ISO_1']: state['properties']['NAME_1'] for state in country_geojson['features']}

        zero_state_counts = pd.DataFrame(columns=['state_name', 'state_code', 'counts'])

        for iso, state_name in geo_states.items():
            if not any(state_counts['state_code'] == iso):  
                new_row = pd.DataFrame({
                    'state_name': [state_name],
                    'state_code': [iso],
                    'counts': [0]
                })
                zero_state_counts = pd.concat([zero_state_counts, new_row], ignore_index=True)

        state_counts = pd.concat([state_counts, zero_state_counts], ignore_index=True)

        state_counts = state_counts.sort_values(by='state_name').reset_index(drop=True)

        state_counts['counts'] = pd.to_numeric(state_counts['counts'], errors='coerce')

        # state_counts.to_csv(f'{selected_country}_state_counts_after_adding_zeroes.csv', index=False)

        fig = px.choropleth(
            state_counts,
            geojson=country_geojson,
            locations='state_code',
            hover_name='state_name',
            hover_data='counts',
            featureidkey='properties.ISO_1',  
            color='counts',
            title=f'Counts for {selected_country}',
            color_continuous_scale=px.colors.sequential.Sunset
        )

        fig.update_layout(
            title=dict(
                text=f"Total isolates from {selected_country}: {len(filtered_data)}<br>Isolates without states' info.: {no_count}",
                xref='paper',
                yref='container',
                y = 0.05,
                x = 0.5,  
            ),
            showlegend=False,
            geo=dict(
                showframe=False,
                showcoastlines=False,
                projection_type='equirectangular'
            ),
            coloraxis_colorbar_title='Counts',

            height=800,
            width=1200,
        )

        fig.update_geos(fitbounds="locations", visible=False)

    return fig

def Assembly_level_bar(df:pd.DataFrame, selected_assembly_levels:list):
    assembly_levels = ['Contig', 'Scaffold', 'Chromosome', 'Complete Genome']
    assembly_level_df = df.groupby('Assembly_level').size().reset_index(name='count')

    assembly_level_df = assembly_level_df.set_index('Assembly_level')
    assembly_level_df = assembly_level_df.reindex(assembly_levels, fill_value=0).reset_index()

    highlight_categories = set(selected_assembly_levels)

    assembly_level_df['selected'] = assembly_level_df['Assembly_level'].isin(highlight_categories)

    fig = px.bar(
        assembly_level_df,
        x='Assembly_level',
        y='count',
        text='count',
        color='selected',
        color_discrete_map={True: 'purple', False: 'lightseagreen'},
        template='plotly_white'
    )

    fig.update_traces(marker=dict(opacity=0.4), selector=dict(marker_color='lightseagreen'))
    fig.update_traces(marker=dict(opacity=0.8), selector=dict(marker_color='purple'))

    fig.update_layout(
        showlegend=False,
        xaxis=dict(
            title=None,
            tickfont=dict(size=11, color='black', family='Arial, sans-serif'),
            tickangle=-45,
            showline= True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        yaxis=dict(
            title='Number of Assemblies',
            titlefont=dict(size=12, color='black', family='Arial, sans-serif'),
            showline=True,
            linecolor='black',
            tickfont=dict(size=11, color='black', family='Arial, sans-serif'),
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        height=400,
        width=400,
    )  

    return fig

def Annotation_bar(df:pd.DataFrame, show_annotations_from:list):

    annotaters = ['GenBank', 'NCBI RefSeq', 'Others', 'No Annotation']

    annotation_df = df.groupby('Annotation_category').size().reset_index(name='count')
    

    annotation_df = annotation_df.set_index('Annotation_category')
    annotation_df = annotation_df.reindex(annotaters, fill_value=0).reset_index()

    
    highlight_categories = set(show_annotations_from)
    annotation_df['selected'] = annotation_df['Annotation_category'].isin(highlight_categories)

    fig = px.bar(
        annotation_df,
        x='Annotation_category',
        y='count',
        text='count',
        color='selected',
        color_discrete_map={True: 'purple', False: 'lightseagreen'},
        template='plotly_white'
    )

    fig.update_traces(marker=dict(opacity=0.4), selector=dict(marker_color='lightseagreen'))
    fig.update_traces(marker=dict(opacity=0.8), selector=dict(marker_color='purple'))


    fig.update_layout(
        showlegend=False,
        xaxis=dict(
            title=None,
            tickfont=dict(size=11, color='black', family='Arial, sans-serif'),
            tickangle=-45,
            showline= True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        yaxis=dict(
            title='Number of Genomes',
            titlefont=dict(size=12, color='black', family='Arial, sans-serif'),
            showline=True,
            linecolor='black',
            tickfont=dict(size=11, color='black', family='Arial, sans-serif'),
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        hovermode='closest',
        hoverlabel=dict(
            bgcolor='white',
            font_size=12,
            font_family='Arial, sans-serif',
            align='left',
        ),
        height=400,
        width=400,
    )

    return fig

def Submission_year_line(df:pd.DataFrame):

    submissions_year_count = df.groupby('Submission_year').size().reset_index(name='count')

    fig = px.line(
    submissions_year_count,
    x='Submission_year',
    y='count',
    template='plotly_white'
    )

    fig.update_traces(mode='lines+markers')

    fig.update_layout(
        xaxis=dict(
            title='Year',
            titlefont=dict(size=12, color='black', family='Arial, sans-serif'),
            tickmode='linear',
            tickfont=dict(size=11, color='black', family='Arial, sans-serif'),
            tick0=1980,
            dtick=1,
            showline=True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        yaxis=dict(
            title='Number of submissions',
            titlefont=dict(size=12, color='black', family='Arial, sans-serif'),
            tickfont=dict(size=11, color='black', family='Arial, sans-serif'),
            showline=True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        plot_bgcolor='rgba(0,0,0,0)',
        height=400,
        width=580,
    )

    return fig

def Sequencing_technologies_scatter(df:pd.DataFrame, selected_sequencing_technologies:list):
    
    plot_df = df.groupby(['Submission_year', 'Categorized_sequencing_technologies']).size().reset_index(name='count')
    
    plot_df.sort_values(by='Categorized_sequencing_technologies', ascending=True)

    plot_df['log_count'] = np.log10(plot_df['count'])
    # plot_df['log_count'] = np.log1p(plot_df['count'])   

    figure = px.scatter(
        plot_df,
        x='Submission_year',
        y='Categorized_sequencing_technologies',
        size='log_count',
        hover_data=['count', 'Categorized_sequencing_technologies'],
        color='Categorized_sequencing_technologies',
        template='plotly_white'
    )

    figure.update_layout(
        xaxis_title="Submission Year",  
        yaxis_title=None,
        showlegend=False,
        hoverlabel=dict(
            bgcolor="white",  
            font_size=12,  
            font_family="Arial",  
        ),
        # margin=dict(t=50, b=50, l=50, r=50),  
        plot_bgcolor='rgba(255,255,255,0)',
        height=700,
        width=930
    )

    figure.update_traces(
        marker=dict(
            opacity=0.7,  
            line=dict(width=1, color='DarkSlateGrey'),  
            sizemode='area',  
            sizeref=2.3 * max(plot_df['log_count'])/(40**2),  
        ),
    )

    return figure

def Coverage_bar(df:pd.DataFrame):

    bin_edges = [0, 50, 100, 200, 300, 400, 500, 700, 900, 1100, 1500, 3000, 5000, float('inf')]
    bin_edges = np.unique(bin_edges)

    df['Coverage_bin'] = pd.cut(df['Coverage_Depth'], bins=bin_edges, right=False)

    full_length = len(df)
    null_length = len(df[pd.isnull(df['Coverage_bin'])])

    coverage_df = df.groupby('Coverage_bin', observed=False).size().reset_index(name='count')

    coverage_df['Coverage_bin'] = coverage_df['Coverage_bin'].astype(str)

    figure = px.bar(
        coverage_df,
        x='Coverage_bin',
        y='count',
        text='count',
        template='plotly_white'
    )


    figure.update_layout(
        showlegend=False,
        title=dict(
                text=f"Total isolates: {full_length}<br>Isolates without coverage data: {null_length}",
                xref='paper',
                yref='container',
                y = 0.85,
                x = 1,  
                font=dict(size=14, color='#202A44'),
        ),
        xaxis=dict(
            title='Coverage Depth (X)',
            titlefont=dict(size=12, color='black', family='Arial, sans-serif'),
            tickfont=dict(size=11, color='black', family='Arial, sans-serif'),
            tickangle=-45,
            showline= True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        yaxis=dict(
            title='Number of Genomes',
            titlefont=dict(size=12, color='black', family='Arial, sans-serif'),
            showline=True,
            linecolor='black',
            tickfont=dict(size=11, color='black', family='Arial, sans-serif'),
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        hovermode='closest',
        hoverlabel=dict(
            bgcolor='white',
            font_size=12,
            font_family='Arial, sans-serif',
            align='left',
        ),
        height=500,
        width=750,
    )

    return figure

def ANI_scatter(df:pd.DataFrame):

    figure = px.scatter(
        df,
        x='ANI_best_match_score',
        y='ANI_best_matched_assembly\'s_coverage',
        symbol='Categorized_sequencing_technologies',
        color='Categorized_sequencing_technologies',
        color_discrete_sequence=px.colors.qualitative.Vivid,
        template='plotly_white'  
    )

    figure.update_layout(
        xaxis=dict(
            title='% Identity',
            titlefont=dict(size=14, color='black', family='Arial, sans-serif'),
            tickmode='linear',
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            tick0=70,
            dtick=5,
            showline=True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        yaxis=dict(
            title='% Coverage',
            titlefont=dict(size=14, color='black', family='Arial, sans-serif'),
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            showline=True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        height=520,
        width=780,
        legend=dict(
            title='Sequencing Technology',
            title_font=dict(size=14, family='Arial, sans-serif'),
            font=dict(size=12, family='Arial, sans-serif'),
            traceorder='normal',
            bgcolor='rgba(255, 255, 255, 0.7)',  
            bordercolor='black',
            borderwidth=1
        )
    )

    return figure

def N50L50_scatter(df:pd.DataFrame):

    figure = px.scatter(
        df,
        x='Contig_N50',
        y='Contig_L50',
        symbol='Categorized_sequencing_technologies',
        color='Categorized_sequencing_technologies',
        color_discrete_sequence=px.colors.qualitative.Vivid,
        template='plotly_white'
    )

    figure.update_layout(
        xaxis=dict(
            title='Scaffold N50 (log scale)',
            titlefont=dict(size=14, color='black', family='Arial, sans-serif'),
            tickmode='linear',
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            showline=True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
            type='log'
            
        ),
        yaxis=dict(
            title='Scaffold L50',
            titlefont=dict(size=14, color='black', family='Arial, sans-serif'),
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            showline=True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
            
        ),
        height=520,
        width=780,
        legend=dict(
            title='Sequencing Technology',
            title_font=dict(size=14, family='Arial, sans-serif'),
            font=dict(size=12, family='Arial, sans-serif'),
            traceorder='normal',
            bgcolor='rgba(255, 255, 255, 0.7)',
            bordercolor='black',
            borderwidth=1
        )
    )

    return figure

def Total_genes_hist(df:pd.DataFrame):

    full_length = len(df)
    null_length = len(df[pd.isnull(df['Total_genes'])])

    figure = px.histogram(
        df,
        x='Total_genes',
        marginal='violin',
        template='plotly_white',
    )

    figure.update_layout(
        showlegend=False,
        title=dict(
                text=f"Total isolates: {full_length}<br>Isolates without 'total gene' data: {null_length}",
                xref='paper',
                yref='container',
                y = 0.90,
                x = 0.5,  
                font=dict(size=10, color='#202A44'),
        ),
        xaxis=dict(
            title='Total Genes',
            titlefont=dict(size=13, color='black', family='Arial, sans-serif'),
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            tickangle=-45,
            showline= True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        yaxis=dict(
            title='Count',
            titlefont=dict(size=13, color='black', family='Arial, sans-serif'),
            showline=True,
            linecolor='black',
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        hovermode='closest',
        hoverlabel=dict(
            bgcolor='white',
            font_size=12,
            font_family='Arial, sans-serif',
            align='left',
        ),
        height=400,
        width=400,
    )

    return figure

def CDSs_hist(df:pd.DataFrame):
    full_length = len(df)
    null_length = len(df[pd.isnull(df['Protein-coding_genes'])])

    figure = px.histogram(
        df,
        x='Protein-coding_genes',
        marginal='violin',
        template='plotly_white',
    )

    figure.update_layout(
        showlegend=False,
        title=dict(
                text=f"Total isolates: {full_length}<br>Isolates without 'CDSs' data: {null_length}",
                xref='paper',
                yref='container',
                y = 0.90,
                x = 0.5,  
                font=dict(size=10, color='#202A44'),
        ),
        xaxis=dict(
            title='Protein Coding genes',
            titlefont=dict(size=13, color='black', family='Arial, sans-serif'),
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            tickangle=-45,
            showline= True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        yaxis=dict(
            title='Count',
            titlefont=dict(size=13, color='black', family='Arial, sans-serif'),
            showline=True,
            linecolor='black',
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        hovermode='closest',
        hoverlabel=dict(
            bgcolor='white',
            font_size=12,
            font_family='Arial, sans-serif',
            align='left',
        ),
        height=400,
        width=400,
    )

    return figure

def Pseudogenes_hist(df:pd.DataFrame):
    full_length = len(df)
    null_length = len(df[pd.isnull(df['Pseudogenes'])])

    figure = px.histogram(
        df,
        x='Pseudogenes',
        marginal='violin',
        template='plotly_white',
    )

    figure.update_layout(
        showlegend=False,
        title=dict(
                text=f"Total isolates: {full_length}<br>Isolates without 'Pseudogene' data: {null_length}",
                xref='paper',
                yref='container',
                y = 0.90,
                x = 0.5,  
                font=dict(size=10, color='#202A44'),
        ),
        xaxis=dict(
            title='Pseudogenes',
            titlefont=dict(size=13, color='black', family='Arial, sans-serif'),
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            tickangle=-45,
            showline= True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        yaxis=dict(
            title='Count',
            titlefont=dict(size=13, color='black', family='Arial, sans-serif'),
            showline=True,
            linecolor='black',
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        hovermode='closest',
        hoverlabel=dict(
            bgcolor='white',
            font_size=12,
            font_family='Arial, sans-serif',
            align='left',
        ),
        height=400,
        width=400,
    )

    return figure

def Non_coding_hist(df:pd.DataFrame):
    full_length = len(df)
    null_length = len(df[pd.isnull(df['Non-coding_genes'])])

    figure = px.histogram(
        df,
        x='Non-coding_genes',
        marginal='violin',
        template='plotly_white',
    )

    figure.update_layout(
        showlegend=False,
        title=dict(
                text=f"Total isolates: {full_length}<br>Isolates without 'Non-coding gene' data: {null_length}",
                xref='paper',
                yref='container',
                y = 0.90,
                x = 0.5,  
                font=dict(size=10, color='#202A44'),
        ),
        xaxis=dict(
            title='Non-coding genes',
            titlefont=dict(size=13, color='black', family='Arial, sans-serif'),
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            tickangle=-45,
            showline= True,
            linecolor='black',
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        yaxis=dict(
            title='Count',
            titlefont=dict(size=13, color='black', family='Arial, sans-serif'),
            showline=True,
            linecolor='black',
            tickfont=dict(size=12, color='black', family='Arial, sans-serif'),
            ticks='outside',
            tickwidth=2,
            ticklen=6,
        ),
        hovermode='closest',
        hoverlabel=dict(
            bgcolor='white',
            font_size=12,
            font_family='Arial, sans-serif',
            align='left',
        ),
        height=400,
        width=400,
    )

    return figure

def Isolation_source_treemap(df:pd.DataFrame):

    if "Unknown" in df['identified_host'].unique():
        columns_of_interest = ['identified_host', 'source_category', 'source', 'sample']
        path = ['identified_host', 'source_category', 'source', 'sample']
    else:
        columns_of_interest = ['source_category', 'source', 'sample']
        
    df.loc[:, columns_of_interest] = df[columns_of_interest].where(pd.notnull(df[columns_of_interest]), "Unknown")
    
    custom_labels = {
        'identified_host': 'Host',
        'source_category': 'Category',
        'source': 'Subcategory',
        'sample': 'Sample'
    }

    fig = px.treemap(
        df, 
        path=['identified_host', 'source_category', 'source', 'sample'],
        width=1000,  
        height=900,
        color='sample',
        hover_data={'source_category': True, 'identified_host': True, 'source': True},
        labels=custom_labels
    )

    return fig
