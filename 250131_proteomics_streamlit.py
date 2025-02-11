import pandas as pd
import streamlit as st

# Load the original protein sequences and TOOL-A/TOOL-B results
uploaded_fasta = st.file_uploader("Upload Original Fasta Data CSV", type=["csv"])
uploaded_tool_a = st.file_uploader("Upload TOOL-A Results CSV", type=["csv"])
uploaded_tool_b = st.file_uploader("Upload TOOL-B Results CSV", type=["csv"])

# Read the CSV files
if uploaded_fasta is not None:
    df = pd.read_csv(uploaded_fasta)
    st.write("Original Fasta Data Preview", df.head())

if uploaded_tool_a is not None:
    tool_a_df = pd.read_csv(uploaded_tool_a)
    st.write("TOOL-A Data Preview", tool_a_df.head())

if uploaded_tool_b is not None:
    tool_b_df = pd.read_csv(uploaded_tool_b)
    st.write("TOOL-B Data Preview", tool_b_df.head())

# Merge the original fasta data with TOOL-A and TOOL-B results based on ProteinName
merged_tool_a = pd.merge(df[['ProteinName', 'Sequence']], tool_a_df[['ProteinName', 'Peptides_A']], on='ProteinName', how='left')
merged_tool_b = pd.merge(df[['ProteinName', 'Sequence']], tool_b_df[['ProteinName', 'Peptides_B']], on='ProteinName', how='left')

# Function to calculate the coverage of a peptide in the protein sequence
def calculate_coverage(protein_sequence, peptide_sequence):
    peptide_length = len(peptide_sequence)
    protein_length = len(protein_sequence)
    
    # Count how many times the peptide appears in the protein sequence
    count = protein_sequence.count(peptide_sequence)
    coverage = count * peptide_length / protein_length
    return coverage

# Initialize lists to store coverage values
coverage_a = []
coverage_b = []

# Iterate over TOOL-A results to calculate coverage
for idx, row in merged_tool_a.iterrows():
    protein_name = row['ProteinName']
    peptide_a = row['Peptides_A']
    full_sequence = row['Sequence']
    
    # Calculate coverage for each peptide found by TOOL-A
    if pd.notna(peptide_a):  # Avoid NaN values
        coverage_a.append(calculate_coverage(full_sequence, peptide_a))
    else:
        coverage_a.append(0)

# Iterate over TOOL-B results to calculate coverage
for idx, row in merged_tool_b.iterrows():
    protein_name = row['ProteinName']
    peptide_b = row['Peptides_B']
    full_sequence = row['Sequence']
    
    # Calculate coverage for each peptide found by TOOL-B
    if pd.notna(peptide_b):  # Avoid NaN values
        coverage_b.append(calculate_coverage(full_sequence, peptide_b))
    else:
        coverage_b.append(0)

# Add coverage columns to merged dataframes
merged_tool_a['Coverage_A'] = coverage_a
merged_tool_b['Coverage_B'] = coverage_b

# Aggregate coverage per ProteinName (average coverage for each protein)
final_aggregated_df = pd.merge(merged_tool_a[['ProteinName', 'Coverage_A']], merged_tool_b[['ProteinName', 'Coverage_B']], on='ProteinName', how='left')

# Calculate average coverage for each tool
final_aggregated_df['Average_Coverage_A'] = final_aggregated_df['Coverage_A']
final_aggregated_df['Average_Coverage_B'] = final_aggregated_df['Coverage_B']

# Final visualization
st.write("Final Aggregated Coverage Data", final_aggregated_df.head())
