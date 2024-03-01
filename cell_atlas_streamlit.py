import pyreadr
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()
from rpy2.robjects.conversion import localconverter
import plotly.graph_objects as go
import numpy as np
import streamlit as st
from pycirclize import Circos
import matplotlib.pyplot as plt



cell_types = ["Adipocytes", "Endothelial Cells", "Erythroblasts", "Fibroblasts",
              "Hematopoietic Cells", "Pericytes", "Schwann Cells",
              "Smooth Muscle Cells", "Vascular Smooth Muscle Cells"]

filepath = 'cell_cell_interactions/cell_atlas_cci_pvat(bat not included)_heatmap_matrices_included.rds' 
readRDS = ro.r['readRDS']
cell_interactions = readRDS(filepath)

st.title('Ligand Receptor Interactions')

source_cell_type = st.selectbox(
    'Select the Source Cell Type',
    cell_types,
    index=0 
)

target_cell_type = st.selectbox(
    'Select the Target Cell Type',
    cell_types,
    index=1  
)

if source_cell_type == target_cell_type:
    
    st.error('Error: The source and target cell type must be different. Please select a different target cell type.')
    
else:
    cell_pair_string = f'{source_cell_type.replace(" ", "")}2{target_cell_type.replace(" ", "")}'
    cell_one2cell_two = cell_interactions.rx2(cell_pair_string)         
    ligand_activities = pandas2ri.rpy2py(cell_one2cell_two.rx2('ligand_activities'))
    ligand_receptor = pandas2ri.rpy2py(cell_one2cell_two.rx2('p_ligand_receptor_network'))
    ligand_target = pandas2ri.rpy2py(cell_one2cell_two.rx2('p_ligand_target_network'))
    filtered_columns = [col for col in ligand_receptor.columns if 'Col' in col]
    filtered_ligand_receptor = ligand_receptor[filtered_columns]
    filtered_ligand_receptor = filtered_ligand_receptor[filtered_ligand_receptor.index.str.contains('Itg')]
    filtered_ligand_receptor = filtered_ligand_receptor.loc[~(filtered_ligand_receptor == 0).all(axis=1)]







st.write("""
         This heatmap visualizes the interactions between ligands acting as source and receptors as target, 
         highlighting the prior strength of each interaction. The intensity of the color 
         indicates the level of interaction, with darker shades representing stronger interactions.
         """)


max_val = ligand_receptor.values.max()
min_val = ligand_receptor[ligand_receptor > 0].min().min()
min_color_intensity = 10 

colorscale = [
    [0, 'white'],  
    [min_color_intensity * min_val / max_val, 'rgb(255,182,193)'], 
    [1, 'rgb(255,105,180)'] 
]


hover_text = [["Ligand: {}<br>Receptor: {}<br>Interaction: {}".format(ligand, receptor, ligand_receptor.loc[receptor, ligand])
               for ligand in ligand_receptor.columns] for receptor in ligand_receptor.index]


fig = go.Figure(data=go.Heatmap(
    z=ligand_receptor.values,  
    x=ligand_receptor.columns,
    y=ligand_receptor.index, 
    text=hover_text,
    hoverinfo='text',
    colorbar=dict(title='Prior Regulatory Potential'),  
    colorscale=colorscale,  
    showscale=True 
))



fig.update_layout(
    #title='Ligand-Receptor Interaction Heatmap',
    xaxis_title='Ligands',
    yaxis_title='Receptors',
    xaxis=dict(
        tickangle=-45,
        side='top',  
        tickfont=dict(size=10),
    ),
    yaxis=dict(
        tickfont=dict(size=10) 
    ),
    plot_bgcolor='white',
    width=1100,
    height=800,
    margin=dict(l=100, r=100, t=100, b=100),  
    paper_bgcolor='white', 
    coloraxis_colorbar=dict(
        thickness=20, 
        len=0.75,  
        yanchor='middle', 
        y=0.5 
    ),
)


st.plotly_chart(fig, use_container_width=True)


st.header('Interactions between Collagens and Integrins')

st.write("From the above heatmap, the below plot has been filtered to show interactions specifically between collagens and integrins.")


circos = Circos.initialize_from_matrix(
    filtered_ligand_receptor.T,
    space=5,
    cmap="tab10",
    label_kws=dict(size=12),
    link_kws=dict(ec="black", lw=0.5, direction=1),
)

fig_circos = circos.plotfig()
st.pyplot(fig_circos)