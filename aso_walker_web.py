import streamlit as st
import pandas as pd
import io
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

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
    # Wallace Rule: 2(A+T) + 4(G+C)
    tm = 2 * (a + t + u) + 4 * (g + c)
    return round(gc_cont, 1), tm

def predict_secondary_structure(sequence):
    """
    Predicts Dot-Bracket using a Thermodynamic-Weighted Nussinov Algorithm.
    Improved to match ViennaRNA behavior (min loop size 3).
    """
    n = len(sequence)
    # Energy weights reflecting stability (Higher = Stronger)
    def energy(b1, b2):
        pair = sorted([b1, b2])
        if pair == ['C', 'G']: return 3.0  # GC is strongest
        if pair == ['A', 'U'] or pair == ['A', 'T']: return 2.0
        if pair == ['G', 'U'] or pair == ['G', 'T']: return 1.0  # Wobble pair
        return 0

    dp = np.zeros((n, n))
    # Fill DP table (Ensuring minimum loop size of 3)
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
            
    # Backtrack to build dot-bracket string
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

def draw_force_directed_structure(sequence, dot_bracket):
    """Generates a Force-Directed 2D plot mimicking professional RNA software."""
    n = len(sequence)
    G = nx.Graph()
    hb_edges = []
    stack = []
    
    # Add backbone and hydrogen bond edges
    for i in range(n - 1):
        G.add_edge(i, i + 1, weight=5)
    
    for i, char in enumerate(dot_bracket):
        if char == "(":
            stack.append(i)
        elif char == ")":
            j = stack.pop()
            G.add_edge(i, j, weight=1)
            hb_edges.append((i, j))
            
    # Physics simulation (Fallback to spring if SciPy is missing)
    try:
        pos = nx.kamada_kawai_layout(G)
    except:
        pos = nx.spring_layout(G, k=1/np.sqrt(n), iterations=50)
    
    fig, ax = plt.subplots(figsize=(10, 10))
    nx.draw_networkx_edges(G, pos, edgelist=[(i, i+1) for i in range(n-1)], 
                           edge_color='#D3D3D3', width=1.5, ax=ax)
    nx.draw_networkx_edges(G, pos, edgelist=hb_edges, 
                           edge_color='#29b09d', width=2.5, ax=ax)
    
    for i, base in enumerate(sequence):
        color = '#0068c9' if base in 'GC' else '#ff4b4b'
        ax.scatter(pos[i][0], pos[i][1], color=color, s=200, zorder=5, edgecolors='white')
        ax.text(pos[i][0], pos[i][1], base, fontsize=8, ha='center', va='center', 
                color='white', weight='bold', zorder=6)
    ax.set_axis_off()
    return fig

# --- UI SECTION ---

st.markdown("<h1 style='text-align: center;'>🧬 ASO Walker Pro</h1>", unsafe_allow_html=True)

# Main Input
raw_seq = st.text_area("Target Sequence (RNA/DNA)", placeholder="Paste sequence here...", height=120)
clean_seq = "".join(raw_seq.upper().split())

# Settings Columns
col_name, col_size, col_step = st.columns(3)
with col_name:
    seq_name = st.text_input("Project Name", value="Task_1")
with col_size:
    aso_size = st.number_input("ASO Size (bp)", min_value=1, value=20)
with col_step:
    step_size = st.slider("Step Size (Stride)", 1, 10, 1)

if st.button("Generate Full Structural Analysis", type="primary", use_container_width=True):
    if not clean_seq:
        st.error("Please enter a sequence.")
    else:
        with st.spinner("Calculating thermodynamic folding..."):
            dot_bracket = predict_secondary_structure(clean_seq)
        
        # 1. Synchronized Viewer
        st.subheader("Interactive Alignment Map")
        ruler = "".join([str(i%10) if i%10==0 else ("1" if i==1 else " ") for i in range(1, len(clean_seq)+1)])
        st.markdown(f"""
        <div style="overflow-x: auto; white-space: pre; font-family: 'Courier New', monospace; 
                    background-color: #f0f2f6; padding: 20px; border-radius: 10px; line-height: 1.6;">
<span style="color:#ff4b4b; font-weight:bold;">NUM:</span> {ruler}
<span style="color:#0068c9; font-weight:bold;">SEQ:</span> {clean_seq}
<span style="color:#29b09d; font-weight:bold;">STR:</span> {dot_bracket}
        </div>""", unsafe_allow_html=True)

        # 2. Folded Topology and Results Table
        col_img, col_data = st.columns([1.2, 1])
        
        with col_img:
            st.subheader("Folded Topology")
            fig = draw_force_directed_structure(clean_seq, dot_bracket)
            st.pyplot(fig)
            st.caption("● Blue: Strong Pairing (GC) | ● Red: Weak Pairing (AU) | Teal: Hydrogen Bonds")

        # 3. Process Table
        results = []
        for i in range(0, len(clean_seq) - aso_size + 1, step_size):
            window_struct = dot_bracket[i : i + aso_size]
            aso_seq = get_reverse_complement(clean_seq[i : i + aso_size])
            gc, tm = calculate_metrics(aso_seq)
            results.append({
                "Name": f"{seq_name}_{len(results) + 1}",
                "Region": f"{i+1}-{i+aso_size}",
                "ASO Sequence (5'-3')": aso_seq,
                "GC%": gc, "Tm (°C)": tm,
                "Accessibility%": round((window_struct.count(".") / aso_size) * 100, 1)
            })
        
        df = pd.DataFrame(results)
        
        with col_data:
            st.subheader("ASO Selection Table")
            st.dataframe(
                df.style.background_gradient(subset=['Accessibility%'], cmap='Greens')
                .format({"GC%": "{:.1f}", "Accessibility%": "{:.1f}"}), 
                use_container_width=True, height=450
            )
            
            csv_buffer = io.StringIO()
            df.to_csv(csv_buffer, index=False)
            st.download_button("💾 Download Analysis CSV", data=csv_buffer.getvalue(), 
                               file_name=f"{seq_name}_Analysis.csv", use_container_width=True)

with st.sidebar:
    st.header("How it works")
    st.info("This app uses a weighted Nussinov algorithm with a 3-base minimum loop constraint to simulate RNA folding.")
    st.divider()
    st.markdown("**Developed for ASO Research**")
