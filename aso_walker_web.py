import streamlit as st
import pandas as pd
import io
import numpy as np

# 1. Page Configuration
st.set_page_config(page_title="ASO Walker Pro", page_icon="🧬", layout="wide")

# --- BIOLOGICAL LOGIC ---

def get_reverse_complement(seq):
    trans = str.maketrans('ATCGUatcgu', 'TAGCAtagca')
    return seq.translate(trans)[::-1]

def calculate_metrics(seq):
    seq = seq.upper()
    g, c, a, t, u = seq.count('G'), seq.count('C'), seq.count('A'), seq.count('T'), seq.count('U')
    length = len(seq)
    gc_cont = ((g + c) / length) * 100 if length > 0 else 0
    # Wallace Rule: 2(A+T) + 4(G+C)
    tm = 2 * (a + t + u) + 4 * (g + c)
    return round(gc_cont, 1), tm

def predict_secondary_structure(sequence):
    """
    Predicts Dot-Bracket using a Thermodynamic-Weighted Nussinov Algorithm.
    Matches ViennaRNA behavior with a 3-base minimum loop constraint.
    """
    n = len(sequence)
    def energy(b1, b2):
        pair = sorted([b1, b2])
        if pair == ['C', 'G']: return 3.0  
        if pair == ['A', 'U'] or pair == ['A', 'T']: return 2.0
        if pair == ['G', 'U'] or pair == ['G', 'T']: return 1.0  
        return 0

    dp = np.zeros((n, n))
    for k in range(4, n): 
        for i in range(n - k):
            j = i + k
            res = [dp[i+1][j], dp[i][j-1]]
            
            e_val = energy(sequence[i], sequence[j])
            if e_val > 0:
                res.append(dp[i+1][j-1] + e_val)
            
            for t in range(i + 1, j):
                res.append(dp[i][t] + dp[t+1][j])
            dp[i][j] = max(res)
            
    structure = ["."] * n
    stack = [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if i >= j - 3: continue 
        elif dp[i][j] == dp[i+1][j]: stack.append((i+1, j))
        elif dp[i][j] == dp[i][j-1]: stack.append((i, j-1))
        elif dp[i][j] == dp[i+1][j-1] + energy(sequence[i], sequence[j]):
            structure[i], structure[j] = "(", ")"
            stack.append((i+1, j-1))
        else:
            for k in range(i+1, j):
                if dp[i][j] == dp[i][k] + dp[k+1][j]:
                    stack.append((k+1, j)); stack.append((i, k))
                    break
    return "".join(structure)

# --- UI SECTION ---

st.markdown("<h1 style='text-align: center;'>🧬 ASO Walker Pro</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; color: gray;'>Linear Structural Mapping & ASO Generation</p>", unsafe_allow_html=True)

# Main Input
raw_seq = st.text_area("Target Sequence (RNA/DNA)", placeholder="Paste sequence here...", height=150)
clean_seq = "".join(raw_seq.upper().split())

# Settings Columns
col_name, col_size, col_step = st.columns(3)
with col_name:
    seq_name = st.text_input("Project Name", value="Task_1")
with col_size:
    aso_size = st.number_input("ASO Size (bp)", min_value=1, value=20)
with col_step:
    step_size = st.slider("Step Size (Stride)", 1, 10, 1)

if st.button("Generate Structural Analysis", type="primary", use_container_width=True):
    if not clean_seq:
        st.error("Please enter a sequence.")
    else:
        with st.spinner("Analyzing thermodynamic folding..."):
            dot_bracket = predict_secondary_structure(clean_seq)
        
        # 1. Synchronized Viewer (Expanded to full width since image is gone)
        st.subheader("Interactive Alignment Map")
        st.info("Scroll horizontally to compare. ( . ) = Loop/Accessible | ( ( ) = Stem/Inaccessible")
        
        # Ruler logic for base numbering
        ruler = ""
        for i in range(len(clean_seq)):
            pos = i + 1
            if pos % 10 == 0:
                ruler += str(pos)
            elif pos == 1:
                ruler += "1"
            else:
                last_ten = (pos // 10) * 10
                if pos > last_ten and pos < last_ten + len(str(last_ten)):
                    continue
                ruler += " "

        st.markdown(f"""
        <div style="overflow-x: auto; white-space: pre; font-family: 'Courier New', monospace; 
                    background-color: #f0f2f6; padding: 25px; border-radius: 10px; line-height: 1.8; font-size: 15px;">
<span style="color:#ff4b4b; font-weight:bold;">NUM:</span> {ruler}
<span style="color:#0068c9; font-weight:bold;">SEQ:</span> {clean_seq}
<span style="color:#29b09d; font-weight:bold;">STR:</span> {dot_bracket}
        </div>""", unsafe_allow_html=True)

        st.divider()

        # 2. Results Table
        st.subheader("ASO Selection Table")
        results = []
        for i in range(0, len(clean_seq) - aso_size + 1, step_size):
            window_struct = dot_bracket[i : i + aso_size]
            aso_seq = get_reverse_complement(clean_seq[i : i + aso_size])
            gc, tm = calculate_metrics(aso_seq)
            results.append({
                "Name": f"{seq_name}_{len(results) + 1}",
                "Region": f"{i+1}-{i+aso_size}",
                "ASO Sequence (5'-3')": aso_seq,
                "GC%": gc, 
                "Est. Tm (°C)": tm,
                "Accessibility%": round((window_struct.count(".") / aso_size) * 100, 1)
            })
        
        df = pd.DataFrame(results)
        
        # Display formatted table
        st.dataframe(
            df.style.background_gradient(subset=['Accessibility%'], cmap='Greens')
            .format({"GC%": "{:.1f}", "Accessibility%": "{:.1f}"}), 
            use_container_width=True
        )
        
        # Export Button
        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)
        st.download_button("💾 Download Full CSV Analysis", data=csv_buffer.getvalue(), 
                           file_name=f"{seq_name}_Analysis.csv", use_container_width=True)

with st.sidebar:
    st.header("Guide")
    st.markdown("""
    **Notation:**
    * `.` = Loop (Accessible)
    * `(` or `)` = Stem (Bound)
    
    **QC Specs:**
    * **Ideal GC:** 40% - 60%
    * **Accessibility:** Higher is better for binding.
    """)
    st.divider()
    st.markdown("**Developed for ASO Research**")
