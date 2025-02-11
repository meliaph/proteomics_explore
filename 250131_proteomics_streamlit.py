import streamlit as st
import pandas as pd

# Upload CSV files
st.title('Proteomics Data Visualization')

st.subheader("Upload the CSV files")

uploaded_tool_a = st.file_uploader("Upload TOOL-A Results CSV", type=["csv"])
uploaded_tool_b = st.file_uploader("Upload TOOL-B Results CSV", type=["csv"])
uploaded_aggregated = st.file_uploader("Upload Final Aggregated Results CSV", type=["csv"])

if uploaded_tool_a is not None:
    tool_a_df = pd.read_csv(uploaded_tool_a)
    st.write("TOOL-A Data Preview", tool_a_df.head())

if uploaded_tool_b is not None:
    tool_b_df = pd.read_csv(uploaded_tool_b)
    st.write("TOOL-B Data Preview", tool_b_df.head())

if uploaded_aggregated is not None:
    final_aggregated_df = pd.read_csv(uploaded_aggregated)
    st.write("Final Aggregated Data Preview", final_aggregated_df.head())

# Function to highlight the peptide sequence in original sequence
def highlight_peptides(original_seq, peptides, color="yellow"):
    highlighted_seq = original_seq
    for peptide in peptides:
        highlighted_seq = highlighted_seq.replace(peptide, f'<mark style="background-color: {color}">{peptide}</mark>')
    return highlighted_seq

# Display the protein sequence with highlighted peptides
if uploaded_tool_a is not None and uploaded_tool_b is not None:
    protein_name = st.text_input("Enter ProteinName to visualize", "1433E_HUMAN")
    
    # Filter the dataframes for the specific ProteinName
    tool_a_peptides = tool_a_df[tool_a_df['ProteinName'] == protein_name]['Peptides_A'].values
    tool_b_peptides = tool_b_df[tool_b_df['ProteinName'] == protein_name]['Peptides_B'].values
    full_sequence = st.text_input("Enter full sequence (if not preloaded):", "")

    # Visualize and highlight peptides for TOOL-A and TOOL-B
    if full_sequence:
        # Highlighting peptides from TOOL-A
        tool_a_highlighted = highlight_peptides(full_sequence, tool_a_peptides, color="yellow")
        st.markdown(f"**TOOL-A Peptides** for {protein_name}:")
        st.markdown(f"<div style='white-space: pre-wrap;'>{tool_a_highlighted}</div>", unsafe_allow_html=True)

        # Highlighting peptides from TOOL-B
        tool_b_highlighted = highlight_peptides(full_sequence, tool_b_peptides, color="lightblue")
        st.markdown(f"**TOOL-B Peptides** for {protein_name}:")
        st.markdown(f"<div style='white-space: pre-wrap;'>{tool_b_highlighted}</div>", unsafe_allow_html=True)

# Visualize aggregated coverage for TOOL-A and TOOL-B
if uploaded_aggregated is not None:
    st.subheader("Coverage Comparison")

    # Plot average coverage for TOOL-A and TOOL-B
    tool_a_avg_coverage = final_aggregated_df['Average_Coverage_A'].mean()
    tool_b_avg_coverage = final_aggregated_df['Average_Coverage_B'].mean()

    st.write(f"**Average Coverage for TOOL-A:** {tool_a_avg_coverage:.5f}")
    st.write(f"**Average Coverage for TOOL-B:** {tool_b_avg_coverage:.5f}")

    # Plot total coverage
    tool_a_total_coverage = final_aggregated_df['Total_Coverage_A'].sum()
    tool_b_total_coverage = final_aggregated_df['Total_Coverage_B'].sum()

    st.write(f"**Total Coverage for TOOL-A:** {tool_a_total_coverage:.5f}")
    st.write(f"**Total Coverage for TOOL-B:** {tool_b_total_coverage:.5f}")
    
    st.bar_chart(final_aggregated_df[['Average_Coverage_A', 'Average_Coverage_B']])

