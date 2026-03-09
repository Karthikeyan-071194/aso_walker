import streamlit as st
import pandas as pd
import io
import numpy as np
import matplotlib.pyplot as plt

# 1. Page Configuration
st.set_page_config(page_title="ASO Walker Pro", page_icon="🧬", layout="wide")

# --- BIOLOGICAL & DRAWING LOGIC ---

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
    """Predicts Dot-Bracket using Nussinov Algorithm."""
    n = len(sequence)
    def energy(b1, b2):
        pair = sorted([b1, b2])
        if pair == ['C', 'G']: return -3
        if pair == ['A', 'U'] or pair == ['A', 'T']: return -2
        if pair == ['G', 'U'] or pair == ['G', 'T']: return -1
        return 0

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
            
    structure = ["."] * n
    stack = [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if i >= j: continue
        elif dp[i][j] == dp[i+1][j]: stack.append((i+1, j))
        elif dp[i][j] == dp[i][j-1]: stack.append((i, j-1))
        elif dp[i][j] == dp[i+1][j-1] + (1 if energy(sequence[i], sequence[j]) < 0 else 0):
            structure[i], structure[j] = "(", ")"
            stack.append((i+1, j-1))
        else:
            for k in range(i+1, j):
                if dp[i][j] == dp[i][k] + dp[k+1][j]:
                    stack.append((k+1, j)); stack.append((i, k))
                    break
    return "".join(structure)

def draw_structure(sequence, dot_bracket):
    """Generates a 2D plot of the folded RNA."""
    n = len(sequence)
    # Simple circular layout for visualization
    angles = np.linspace(0, 2*np.pi, n, endpoint=False)
    x = np.cos(angles)
    y = np.sin(angles)
    
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(x, y, color='lightgrey', lw=1, zorder=1) # Backbone
    
    # Draw base pairs
    stack = []
    for i, char in enumerate(dot_bracket):
        if char == "(":
            stack.append(i)
        elif char == ")":
            j = stack.pop()
            ax.plot([x[i], x[j]], [y[i], y[j]], color='#29b09d', lw=2, alpha=0.6, zorder=2)
            
    # Draw nucleotides
    for i, base in enumerate(sequence):
        color = '#0068c9' if base in 'GC' else '#ff4b4b'
        ax.scatter(x[i], y[i], color=color, s=100, zorder=3)
        ax.text(x[i], y[i], base, fontsize=8, ha='center', va='center', color='white', weight='bold', zorder=4)

    ax.set_axis_off()
    ax.set_title("Predicted 2D Folding Map", fontsize=10, family='monospace')
    return fig

# --- UI SECTION ---

st.markdown("<h1 style='text-align: center;'>🧬 ASO Walker Pro</h1>", unsafe_allow_html=True)

raw_seq = st.text_area("Target Sequence", placeholder="Paste sequence here...", height=100)
clean_seq = "".join(raw_seq.upper().split())

col_p, col_s, col_st = st.columns(3)
with col_p:
    seq_name = st.text_input("Project Name", value="Task_1")
with col_s:
    aso_size = st.number_input("ASO Size (bp)", min_value=1, value=20)
with col_st:
    step_size = st.slider("Step Size", 1, 10, 1)

if st.button("Generate Folded Image & ASOs", type="primary", use_container_width=True):
    if not clean_seq:
        st.error("Please enter a sequence.")
    else:
        # 1. Prediction
        with st.spinner("Calculating folding energy..."):
            dot_bracket = predict_secondary_structure(clean_seq)
        
        # 2. Visual Layout
        col_map, col_img = st.columns([2, 1])
        
        with col_map:
            st.subheader("Alignment Map")
            ruler = "".join([str(i%10) if i%10==0 else ("1" if i==1 else " ") for i in range(1, len(clean_seq)+1)])
            st.markdown(f"""
            <div style="overflow-x: auto; white-space: pre; font-family: 'Courier New', monospace; 
                        background-color: #f0f2f6; padding: 15px; border-radius: 5px; font-size: 14px;">
<span style="color:red">NUM:</span> {ruler}
<span style="color:blue">SEQ:</span> {clean_seq}
<span style="color:teal">STR:</span> {dot_bracket}
            </div>""", unsafe_allow_html=True)

        with col_img:
            st.subheader("2D Structure")
            fig = draw_structure(clean_seq, dot_bracket)
            st.pyplot(fig)

        # 3. Table
        results = []
        for i in range(0, len(clean_seq) - aso_size + 1, step_size):
            window_struct = dot_bracket[i : i + aso_size]
            aso_seq = get_reverse_complement(clean_seq[i : i + aso_size])
            gc, tm = calculate_metrics(aso_seq)
            results.append({
                "Name": f"{seq_name}_{len(results) + 1}",
                "Sequence 5'-3'": aso_seq,
                "GC%": gc, "Tm (°C)": tm,
                "Accessibility%": round((window_struct.count(".") / aso_size) * 100, 1)
            })
        
        df = pd.DataFrame(results)
        st.dataframe(df.style.background_gradient(subset=['Accessibility%'], cmap='Greens'), use_container_width=True)
        
        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)
        st.download_button("💾 Download Full Analysis", data=csv_buffer.getvalue(), file_name=f"{seq_name}_Analysis.csv", use_container_width=True)

with st.sidebar:
    st.info("● Blue: G/C (Strong) | ● Red: A/U (Weak) | Green Lines: Hydrogen Bonds")
    st.write("Developed for ASO Research")
