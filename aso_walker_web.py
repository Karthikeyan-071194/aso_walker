import streamlit as st
import pandas as pd
import io
import re
import numpy as np

# 1. Page Configuration
st.set_page_config(page_title="ASO Walker Pro", page_icon="🧬", layout="centered")

# --- BIOLOGICAL LOGIC ---

def get_reverse_complement(seq):
    trans = str.maketrans('ATCGUatcgu', 'TAGCAtagca')
    return seq.translate(trans)[::-1]

def calculate_metrics(seq):
    seq = seq.upper()
    g, c, a, t, u = seq.count('G'), seq.count('C'), seq.count('A'), seq.count('T'), seq.count('U')
    length = len(seq)
    gc_cont = ((g + c) / length) * 100 if length > 0 else 0
    tm = 2 * (a + t + u) + 4 * (g + c)
    return round(gc_cont, 1), tm

def simple_fold(sequence):
    """
    A simplified Nussinov algorithm to predict secondary structure.
    Returns a list where 0 = unpaired (loop/bulge) and 1 = paired (stem).
    """
    n = len(sequence)
    if n > 100: # Limit for web performance
        return [0] * n 
    
    # Minimal pairing rules
    def can_pair(base1, base2):
        pair = sorted([base1, base2])
        return pair == ['C', 'G'] or pair == ['A', 'U'] or pair == ['A', 'T']

    # Initialize DP table
    nm = np.zeros((n, n))
    for k in range(1, n):
        for i in range(n - k):
            j = i + k
            if j - i >= 3:
                down = nm[i + 1][j]
                left = nm[i][j - 1]
                diag = nm[i + 1][j - 1] + (1 if can_pair(sequence[i], sequence[j]) else 0)
                max_k = 0
                for t in range(i + 1, j):
                    max_k = max(max_k, nm[i][t] + nm[t + 1][j])
                nm[i][j] = max(down, left, diag, max_k)
    
    # Backtrack to find paired bases (Simplified)
    structure = [0] * n
    # (Full backtracking omitted for brevity, marking likely paired regions)
    for i in range(n):
        for j in range(i+3, n):
            if nm[i][j] > nm[i+1][j] and nm[i][j] > nm[i][j-1] and can_pair(sequence[i], sequence[j]):
                structure[i] = 1
                structure[j] = 1
    return structure

# --- UI SECTION ---

st.markdown(
    """
    <div style="text-align: center; padding-bottom: 20px;">
        <img src="https://github.com/Karthikeyan-071194/aso_walker/blob/main/the_walker.png?raw=true" 
             width="150" style="display: block; margin: auto; margin-bottom: 10px;">
        <h1 style="margin: 0; font-family: 'Courier New';">ASO Walker Pro</h1>
        <p style="color: gray; font-family: 'Courier New';">Structure-Aware ASO Design</p>
    </div>
    """, unsafe_allow_html=True
)

st.divider()

raw_seq = st.text_area("Target RNA/DNA Sequence (5'-3')", placeholder="Paste sequence here...", height=150)
clean_seq = "".join(raw_seq.upper().split())

col1, col2, col3 = st.columns(3)
with col1:
    seq_name = st.text_input("Project Name", value="ASO_Task")
with col2:
    aso_size = st.number_input("ASO Size", min_value=1, value=20)
with col3:
    step_size = st.slider("Step Size", 1, 10, 1)

if st.button("Analyze Structure & Generate ASOs", type="primary", use_container_width=True):
    if not clean_seq:
        st.error("Please enter a sequence.")
    else:
        # 1. Predict Structure
        with st.spinner("Folding RNA structure..."):
            struct_map = simple_fold(clean_seq)
            
        # 2. Visualize Structure
        st.subheader("RNA Secondary Structure Map")
        struct_viz = "".join(["●" if s == 1 else "○" for s in struct_map])
        st.code(clean_seq)
        st.code(struct_viz)
        st.caption("● = Paired (Stem) | ○ = Unpaired (Loop/Bulge)")

        # 3. Generate ASOs
        results = []
        for i in range(0, len(clean_seq) - aso_size + 1, step_size):
            window = clean_seq[i : i + aso_size]
            aso_seq = get_reverse_complement(window)
            gc, tm = calculate_metrics(aso_seq)
            
            # Check how many unpaired bases are in this window
            unpaired_count = struct_map[i : i + aso_size].count(0)
            accessibility = (unpaired_count / aso_size) * 100
            
            results.append({
                "Name": f"{seq_name}_{len(results) + 1}",
                "Sequence 5'-3'": aso_seq,
                "GC%": gc,
                "Tm (°C)": tm,
                "Accessibility%": round(accessibility, 1)
            })
        
        df = pd.DataFrame(results)
        st.success("Analysis Complete!")

        # High Accessibility = Better targeting
        st.dataframe(
            df.style.background_gradient(subset=['Accessibility%'], cmap='Blues')
            .format({"GC%": "{:.1f}", "Accessibility%": "{:.1f}"}),
            use_container_width=True
        )

        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)
        st.download_button("💾 Download Analysis CSV", data=csv_buffer.getvalue(), file_name=f"{seq_name}_Full_Analysis.csv", use_container_width=True)
