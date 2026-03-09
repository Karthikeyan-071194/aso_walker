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

def predict_secondary_structure(sequence):
    """
    Predicts RNA secondary structure using a simplified energy-based folding model.
    Returns Dot-Bracket notation: '.' (Loop), '(' or ')' (Stem).
    """
    n = len(sequence)
    # Energy table for pairings
    def energy(b1, b2):
        pair = sorted([b1, b2])
        if pair == ['C', 'G']: return -3
        if pair == ['A', 'U'] or pair == ['A', 'T']: return -2
        if pair == ['G', 'U'] or pair == ['G', 'T']: return -1
        return 0

    # Fill DP table (Nussinov-Jacobson variation)
    dp = np.zeros((n, n))
    for k in range(1, n):
        for i in range(n - k):
            j = i + k
            res = [dp[i+1][j], dp[i][j-1]]
            if energy(sequence[i], sequence[j]) < 0:
                res.append(dp[i+1][j-1] + 1)
            for t in range(i+1, j):
                res.append(dp[i][t] + dp[t+1][j])
            dp[i][j] = max(res)
            
    # Backtrack to build dot-bracket string
    structure = ["."] * n
    stack = [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if i >= j: continue
        elif dp[i][j] == dp[i+1][j]:
            stack.append((i+1, j))
        elif dp[i][j] == dp[i][j-1]:
            stack.append((i, j-1))
        elif dp[i][j] == dp[i+1][j-1] + (1 if energy(sequence[i], sequence[j]) < 0 else 0):
            structure[i] = "("
            structure[j] = ")"
            stack.append((i+1, j-1))
        else:
            for k in range(i+1, j):
                if dp[i][j] == dp[i][k] + dp[k+1][j]:
                    stack.append((k+1, j))
                    stack.append((i, k))
                    break
    return "".join(structure)

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
        # 1. Fold the RNA
        with st.spinner("Folding RNA to find loops and bulges..."):
            dot_bracket = predict_secondary_structure(clean_seq)
            
        # 2. Visual Representation of Folding
        st.subheader("RNA Secondary Structure Prediction")
        
        
        
        st.info("The structure below shows how the RNA 'folds' on itself. Focus your ASO design on the **Dots (.)** which represent accessible loops.")
        
        # Display Dot-Bracket Map
        st.code(f"SEQ: {clean_seq}")
        st.code(f"STR: {dot_bracket}")
        
        # 3. Generate ASOs with Accessibility Analysis
        results = []
        for i in range(0, len(clean_seq) - aso_size + 1, step_size):
            window_seq = clean_seq[i : i + aso_size]
            window_struct = dot_bracket[i : i + aso_size]
            
            aso_seq = get_reverse_complement(window_seq)
            gc, tm = calculate_metrics(aso_seq)
            
            # Accessibility: Percentage of target region that is NOT a stem
            unpaired_bases = window_struct.count(".")
            accessibility = (unpaired_bases / aso_size) * 100
            
            results.append({
                "Name": f"{seq_name}_{len(results) + 1}",
                "Target Region": f"{i+1}-{i+aso_size}",
                "Sequence 5'-3'": aso_seq,
                "GC Content (%)": gc,
                "Est. Tm (°C)": tm,
                "Accessibility (%)": round(accessibility, 1)
            })
        
        df = pd.DataFrame(results)
        st.success("Analysis Complete!")

        # High Accessibility = Darker Blue (Target these!)
        st.dataframe(
            df.style.background_gradient(subset=['Accessibility (%)'], cmap='Blues')
            .format({"GC Content (%)": "{:.1f}", "Accessibility (%)": "{:.1f}"}),
            use_container_width=True
        )

        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)
        st.download_button("💾 Download Full Structural Analysis", data=csv_buffer.getvalue(), file_name=f"{seq_name}_Analysis.csv", use_container_width=True)

with st.sidebar:
    st.header("Structure Guide")
    st.write("**( . ) Dot:** Unpaired loop/bulge (High Binding)")
    st.write("**( ( ) Bracket:** Paired stem (Low Binding)")
    st.divider()
    st.write("Developed for ASO Research")

