import streamlit as st
import pandas as pd
import re
import matplotlib.pyplot as plt
from functools import reduce

# Streamlit app title
st.title("Protein Sequence Visualization with Simplified-Semi-Tryptic Classification")

# File uploaders
main_file = st.file_uploader("Upload Main Data CSV", type=["csv"])
tool_a_file = st.file_uploader("Upload TOOL-A CSV", type=["csv"])
tool_b_file = st.file_uploader("Upload TOOL-B CSV", type=["csv"])

if main_file and tool_a_file and tool_b_file:
    # Read CSV files
    main_df = pd.read_csv(main_file)
    tool_a_df = pd.read_csv(tool_a_file)
    tool_b_df = pd.read_csv(tool_b_file)

    # Display data samples
    st.subheader("Main Data Sample")
    st.dataframe(main_df.head())

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("TOOL-A Sample")
        st.dataframe(tool_a_df.head())
    with col2:
        st.subheader("TOOL-B Sample")
        st.dataframe(tool_b_df.head())

    # Merge data on ProteinName
    merged_df = reduce(lambda left, right: pd.merge(left, right, on="ProteinName", how="left"), 
                       [main_df, tool_a_df, tool_b_df])
    
    # Dropdown for protein selection
    selected_protein = st.selectbox("Select a Protein Name", merged_df["ProteinName"].dropna().unique())

    if selected_protein:
        # Extract sequence and peptides
        sequence = merged_df.loc[merged_df["ProteinName"] == selected_protein, "Sequence"].values[0]
        peptides_a = merged_df.loc[merged_df["ProteinName"] == selected_protein, "Peptides_A"].dropna().unique()
        peptides_b = merged_df.loc[merged_df["ProteinName"] == selected_protein, "Peptides_B"].dropna().unique()

        # Classification function
        def classify_simplified_semi_tryptic(protein_sequence, peptide):
            if pd.isna(protein_sequence) or pd.isna(peptide):
                return False
            
            index = protein_sequence.find(peptide)
            if index == -1:
                return False

            if peptide[-1] in ['K', 'R']:
                if index == 0 or protein_sequence[index - 1] not in ['K', 'R']:
                    return True
            
            if peptide[-1] not in ['K', 'R']:
                if index > 0 and protein_sequence[index - 1] in ['K', 'R']:
                    return True

            return False

        # Toggle buttons for highlighting
        if "highlight_a" not in st.session_state:
            st.session_state["highlight_a"] = False
        if "highlight_b" not in st.session_state:
            st.session_state["highlight_b"] = False

        if st.button("Highlight Peptides_A"):
            st.session_state["highlight_a"] = not st.session_state["highlight_a"]
        if st.button("Highlight Peptides_B"):
            st.session_state["highlight_b"] = not st.session_state["highlight_b"]

        # Highlighting function
        def highlight_sequence(seq, peptides_a, peptides_b):
            for pep in peptides_a:
                if isinstance(pep, str):
                    style = "text-decoration: underline;" if classify_simplified_semi_tryptic(seq, pep) else "background-color:red;"
                    seq = re.sub(f"({pep})", f"<span style='{style}'>{pep}</span>", seq)
            for pep in peptides_b:
                if isinstance(pep, str):
                    style = "text-decoration: underline;" if classify_simplified_semi_tryptic(seq, pep) else "color:blue; font-weight:bold;"
                    seq = re.sub(f"({pep})", f"<span style='{style}'>{pep}</span>", seq)
            return seq

        # Apply highlighting
        highlighted_seq = sequence
        if st.session_state["highlight_a"] or st.session_state["highlight_b"]:
            highlighted_seq = highlight_sequence(sequence, peptides_a if st.session_state["highlight_a"] else [], 
                                                 peptides_b if st.session_state["highlight_b"] else [])

        # Highlight K and R responsible for classification in red
        highlighted_seq = re.sub(r'([KR])', r"<span style='color:red; font-weight:bold;'>\1</span>", highlighted_seq)

        # Display highlighted sequence
        st.subheader("Highlighted Protein Sequence")
        st.markdown(f"""<div style='font-family:monospace; font-size:18px; white-space:pre-wrap; word-wrap:break-word;'>{highlighted_seq}</div>""", unsafe_allow_html=True)
        
        # Create a table displaying Peptides, Coverage, and Sequence Length
        peptide_data = []
        for pep in peptides_a:
            coverage = (len(pep) / len(sequence)) * 100
            peptide_data.append([selected_protein, pep, coverage, len(sequence)])
        for pep in peptides_b:
            coverage = (len(pep) / len(sequence)) * 100
            peptide_data.append([selected_protein, pep, coverage, len(sequence)])
        
        df_peptide_data = pd.DataFrame(peptide_data, columns=["ProteinName", "Peptide", "Coverage", "Sequence Length"])
        
        st.subheader("Peptide Coverage Table")
        st.dataframe(df_peptide_data)
