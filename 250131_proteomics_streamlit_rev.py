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
                return "no_simplified_semi_tryptic"
            
            index = protein_sequence.find(peptide)
            if index == -1:
                return "no_simplified_semi_tryptic"

            if peptide[-1] in ['K', 'R']:
                if index == 0 or protein_sequence[index - 1] not in ['K', 'R']:
                    return "simplified_semi_tryptic"
            
            if peptide[-1] not in ['K', 'R']:
                if index > 0 and protein_sequence[index - 1] in ['K', 'R']:
                    return "simplified_semi_tryptic"

            return "no_simplified_semi_tryptic"

        # Classify peptides
        simplified_semi_tryptic_counts = {"TOOL-A": 0, "TOOL-B": 0}
        
        def classify_and_count(peptides, tool):
            count = 0
            for pep in peptides:
                if classify_simplified_semi_tryptic(sequence, pep) == "simplified_semi_tryptic":
                    count += 1
            simplified_semi_tryptic_counts[tool] = count

        classify_and_count(peptides_a, "TOOL-A")
        classify_and_count(peptides_b, "TOOL-B")

        # Display coverage
        st.subheader("Coverage")
        st.write(f"Coverage TOOL-A: {simplified_semi_tryptic_counts['TOOL-A']:.2f}%")
        st.write(f"Coverage TOOL-B: {simplified_semi_tryptic_counts['TOOL-B']:.2f}%")
        st.write(f"Total Coverage: {simplified_semi_tryptic_counts['TOOL-A'] + simplified_semi_tryptic_counts['TOOL-B']:.2f}%")
        
        # Coverage visualization
        fig, ax = plt.subplots(figsize=(18, 2))
        bars = ax.barh(["TOOL-A", "TOOL-B", "Total"], 
                        [simplified_semi_tryptic_counts["TOOL-A"], 
                         simplified_semi_tryptic_counts["TOOL-B"], 
                         simplified_semi_tryptic_counts["TOOL-A"] + simplified_semi_tryptic_counts["TOOL-B"]], 
                        color=["red", "blue", "green"])
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 1, bar.get_y() + bar.get_height()/2, f'{width:.0f}', va='center')
        ax.set_xlabel("Count")
        ax.set_title("Simplified-Semi-Tryptic Peptide Count")
        st.pyplot(fig)
