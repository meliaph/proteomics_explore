# -*- coding: utf-8 -*-
"""250131_Proteomics_Streamlit.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1NYpQVzoJXtnj0pb90D12dQ7R5MZ2GtRy
"""

import streamlit as st
import pandas as pd

# Load the CSV data from GitHub (or local file for testing)
df = pd.read_csv('your_github_link_or_local_file.csv')

# Streamlit interface
st.title('Protein Peptide Coverage Visualization')

# Dropdown to select ProteinName
protein_name = st.selectbox('Select Protein', df['ProteinName'].unique())

# Filter data for the selected protein
protein_data = df[df['ProteinName'] == protein_name].iloc[0]

# Get the full sequence and peptides
sequence = protein_data['Sequence']
peptides_a = protein_data['Peptides_A'].split(',')  # Assuming Peptides_A is a comma-separated string
peptides_b = protein_data['Peptides_B'].split(',')

# Display the full sequence with highlighted peptides
highlighted_sequence_a = sequence
highlighted_sequence_b = sequence

# Highlight peptides in the sequence (this is a simple example, you could use more advanced regex or string manipulation)
for peptide in peptides_a:
    highlighted_sequence_a = highlighted_sequence_a.replace(peptide, f"<span style='background-color: yellow'>{peptide}</span>")

for peptide in peptides_b:
    highlighted_sequence_b = highlighted_sequence_b.replace(peptide, f"<span style='background-color: lightblue'>{peptide}</span>")

# Display the protein sequence with highlighted peptides
st.markdown(f"**Protein Sequence for {protein_name}:**")
st.markdown(f"<div>{highlighted_sequence_a}</div>", unsafe_allow_html=True)

# Display the coverage for TOOL-A and TOOL-B
coverage_a = protein_data['Coverage_A']
coverage_b = protein_data['Coverage_B']
st.write(f"**Coverage for TOOL-A:** {coverage_a}%")
st.write(f"**Coverage for TOOL-B:** {coverage_b}%")

# Optionally, you can use bar charts to show comparison
st.bar_chart([coverage_a, coverage_b], use_container_width=True, height=200)