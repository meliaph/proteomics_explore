import streamlit as st
import pandas as pd

# Streamlit interface
st.title('Proteomics Data Visualization')

# Upload CSV files
st.subheader("Upload the CSV files")

uploaded_tool_a = st.file_uploader("Upload TOOL-A Results CSV", type=["csv"])
uploaded_tool_b = st.file_uploader("Upload TOOL-B Results CSV", type=["csv"])
uploaded_aggregated = st.file_uploader("Upload Final Aggregated Results CSV", type=["csv"])

# Check if files are uploaded and read them
if uploaded_tool_a is not None:
    tool_a_df = pd.read_csv(uploaded_tool_a)
    st.write("TOOL-A Data Preview", tool_a_df.head())

if uploaded_tool_b is not None:
    tool_b_df = pd.read_csv(uploaded_tool_b)
    st.write("TOOL-B Data Preview", tool_b_df.head())

if uploaded_aggregated is not None:
    final_aggregated_df = pd.read_csv(uploaded_aggregated)
    st.write("Final Aggregated Data Preview", final_aggregated_df.head())

# If the required data is uploaded, proceed with visualization
if uploaded_tool_a is not None and uploaded_tool_b is not None and uploaded_aggregated is not None:
    # Merge dataframes with the final aggregated result
    merged_tool_a = pd.merge(final_aggregated_df[['ProteinName', 'Coverage_A']], tool_a_df[['ProteinName', 'Peptides_A']], on='ProteinName', how='left')
    merged_tool_b = pd.merge(final_aggregated_df[['ProteinName', 'Coverage_B']], tool_b_df[['ProteinName', 'Peptides_B']], on='ProteinName', how='left')

    # Dropdown to select ProteinName
    protein_name = st.selectbox('Select Protein', final_aggregated_df['ProteinName'].unique())

    # Filter data for the selected protein
    protein_data = final_aggregated_df[final_aggregated_df['ProteinName'] == protein_name].iloc[0]

    # Get the full sequence and peptides for TOOL-A and TOOL-B
    sequence = protein_data['Sequence']
    peptides_a = merged_tool_a[merged_tool_a['ProteinName'] == protein_name]['Peptides_A'].values[0].split(',')
    peptides_b = merged_tool_b[merged_tool_b['ProteinName'] == protein_name]['Peptides_B'].values[0].split(',')

    # Display the full protein sequence
    st.markdown(f"**Protein Sequence for {protein_name}:**")
    st.text(sequence)

    # Option to display peptides from TOOL-A, TOOL-B, or both
    show_tool_a = st.checkbox('Show TOOL-A Peptides')
    show_tool_b = st.checkbox('Show TOOL-B Peptides')

    highlighted_sequence = sequence

    # Highlight peptides in the sequence
    if show_tool_a:
        for peptide in peptides_a:
            highlighted_sequence = highlighted_sequence.replace(peptide, f"<span style='background-color: yellow'>{peptide}</span>")

    if show_tool_b:
        for peptide in peptides_b:
            highlighted_sequence = highlighted_sequence.replace(peptide, f"<span style='background-color: lightblue'>{peptide}</span>")

    # Display the protein sequence with highlighted peptides
    st.markdown(f"<div>{highlighted_sequence}</div>", unsafe_allow_html=True)

    # Display the coverage for TOOL-A and TOOL-B
    coverage_a = protein_data['Coverage_A']
    coverage_b = protein_data['Coverage_B']

    st.write(f"**Coverage for TOOL-A:** {coverage_a}%")
    st.write(f"**Coverage for TOOL-B:** {coverage_b}%")

    # Optionally, you can use bar charts to show comparison
    st.bar_chart([coverage_a, coverage_b], use_container_width=True, height=200)
