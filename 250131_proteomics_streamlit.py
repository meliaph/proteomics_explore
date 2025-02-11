import streamlit as st
import pandas as pd
import re
import matplotlib.pyplot as plt
from functools import reduce

# Streamlit app title
st.title("Protein Sequence Visualization")

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
        
        # Initialize state for highlighting
        if "highlight_a" not in st.session_state:
            st.session_state["highlight_a"] = False
        if "highlight_b" not in st.session_state:
            st.session_state["highlight_b"] = False
        
        # Toggle buttons for highlighting
        if st.button("Highlight Peptides_A"):
            st.session_state["highlight_a"] = not st.session_state["highlight_a"]
        if st.button("Highlight Peptides_B"):
            st.session_state["highlight_b"] = not st.session_state["highlight_b"]
        
        # Function to highlight sequence
        def highlight_sequence(seq, peptides_a, peptides_b):
            for pep in peptides_a:
                if isinstance(pep, str):
                    seq = re.sub(f"({pep})", f"<span style='background-color:red;'>{pep}</span>", seq)
            for pep in peptides_b:
                if isinstance(pep, str):
                    seq = re.sub(f"({pep})", f"<span style='color:blue; font-weight:bold;'>{pep}</span>", seq)
            return seq
        
        # Apply highlighting
        highlighted_seq = sequence
        if st.session_state["highlight_a"] or st.session_state["highlight_b"]:
            highlighted_seq = highlight_sequence(sequence, peptides_a if st.session_state["highlight_a"] else [], 
                                                 peptides_b if st.session_state["highlight_b"] else [])
        
        # Coverage calculation (Avoid Overlapping Counts)
        def calculate_coverage(seq, peptides):
            covered_positions = set()
            for pep in peptides:
                if isinstance(pep, str):
                    start_idx = 0
                    while (start_idx := seq.find(pep, start_idx)) != -1:
                        covered_positions.update(range(start_idx, start_idx + len(pep)))
                        start_idx += 1  # Move forward to find other occurrences
            return (len(covered_positions) / len(seq)) * 100 if len(seq) > 0 else 0
        
        coverage_a = calculate_coverage(sequence, peptides_a) if st.session_state["highlight_a"] else 0
        coverage_b = calculate_coverage(sequence, peptides_b) if st.session_state["highlight_b"] else 0
        total_coverage = calculate_coverage(sequence, list(peptides_a) + list(peptides_b))

        # Display sequence with highlighting
        st.subheader("Highlighted Protein Sequence")
        st.markdown(f"""<div style='font-family:monospace; font-size:18px; white-space:pre-wrap; word-wrap:break-word;'>{highlighted_seq}</div>""", unsafe_allow_html=True)
        
         # Display coverage
        st.subheader("Coverage")
        st.write(f"Coverage TOOL-A: {coverage_a:.2f}%")
        st.write(f"Coverage TOOL-B: {coverage_b:.2f}%")
        st.write(f"Total Coverage: {total_coverage:.2f}%")
        
        # Coverage visualization (horizontal bar chart)
        fig, ax = plt.subplots(figsize=(18, 6))
        bars = ax.barh(["TOOL-A", "TOOL-B", "Total"], [coverage_a, coverage_b, total_coverage], color=["red", "blue", "green"])
        
        # Annotate the bars with their values
        for bar in bars:
            width = bar.get_width()
            ax.text(width + 1, bar.get_y() + bar.get_height()/2, f'{width:.2f}%', va='center')
        
        ax.set_xlabel("Coverage (%)")
        ax.set_title("Coverage Comparison")
        ax.set_xlim(0, 100)
        st.pyplot(fig)

             
        

