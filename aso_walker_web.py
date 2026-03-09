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
    tm = 2 * (a + t + u) + 4 * (g + c)
    return round(gc_cont, 1), tm

def predict_secondary_structure(sequence):
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

# --- UI SECTION ---

st.markdown(
    """
    <div style="text-align: center; padding-bottom: 10px;">
        <h1 style="margin: 0; font-family: 'Courier New';">🧬 ASO Walker Pro</h1>
        <p style="color: gray;">Synchronized Structure & Sequence Map</p>
    </div>
    """, unsafe_allow_html=True
)

raw_seq = st.text_area("Target RNA/DNA Sequence", placeholder="Paste sequence here...", height=120)
clean_seq = "".join(raw_seq.upper().split())

col1, col2, col3 = st.columns(3)
with col1:
    seq_name = st.text_input("Project Name", value="ASO_Task")
with col2:
    aso_size = st.number_input("ASO Size (bp)", min_value=1, value=20)
with col3:
    step_size = st.slider("Step Size (Stride)", 1, 10, 1)

if st.button("Generate Synchronized Map", type="primary", use_container_width=True):
    if not clean_seq:
        st.error("Please enter a sequence.")
    else:
        # 1. Fold logic
        with st.spinner("Analyzing folding patterns..."):
            dot_bracket = predict_secondary_structure(clean_seq)
        
        # 2. Build the Synchronized Ruler (Numbers)
        # We create a string where numbers appear every 10 bases
        num_line = ""
        for i in range(1, len(clean_seq) + 1):
            if i == 1 or i % 10 == 0:
                num_line += str(i).ljust(1)
            else:
                # Add a placeholder dot for spacing
                if len(str(i-1)) > 1 and (i-1) % 10 == 0:
                     continue # Handle multi-digit spacing
                num_line += " "
        
        # Custom logic to ensure numbers align with characters
        ruler = ""
        for i in range(len(clean_seq)):
            pos = i + 1
            if pos % 10 == 0:
                ruler += str(pos)
            elif pos == 1:
                ruler += "1"
            else:
                # Only add space if we aren't currently inside a multi-digit number
                last_ten = (pos // 10) * 10
                if pos > last_ten and pos < last_ten + len(str(last_ten)):
                    continue
                ruler += " "

        # --- THE SYNCHRONIZED SCROLL WINDOW ---
        st.subheader("Interactive Alignment Map")
        st.info("Scroll horizontally to compare. ( . ) = Loop/Accessible | ( ( ) = Stem/Inaccessible")
        
        # We use a white-space: pre CSS tag to keep everything monospaced and locked
        st.markdown(
            f"""
            <div style="overflow-x: auto; white-space: pre; font-family: 'Courier New', monospace; 
                        background-color: #f0f2f6; padding: 20px; border-radius: 10px; line-height: 1.5; color: #31333F;">
<span style="color: #ff4b4b; font-weight: bold;">NUM:</span> {ruler}
<span style="color: #0068c9; font-weight: bold;">SEQ:</span> {clean_seq}
<span style="color: #29b09d; font-weight: bold;">STR:</span> {dot_bracket}
            </div>
            """, unsafe_allow_html=True
        )

        # 3. Table Generation
        results = []
        for i in range(0, len(clean_seq) - aso_size + 1, step_size):
            window_struct = dot_bracket[i : i + aso_size]
            aso_seq = get_reverse_complement(clean_seq[i : i + aso_size])
            gc, tm = calculate_metrics(aso_seq)
            accessibility = (window_struct.count(".") / aso_size) * 100
            
            results.append({
                "Name": f"{seq_name}_{len(results) + 1}",
                "Target Region": f"{i+1}-{i+aso_size}",
                "Sequence 5'-3'": aso_seq,
                "GC%": gc,
                "Est. Tm (°C)": tm,
                "Accessibility%": round(accessibility, 1)
            })
        
        df = pd.DataFrame(results)
        st.dataframe(
            df.style.background_gradient(subset=['Accessibility%'], cmap='Greens')
            .format({"GC%": "{:.1f}", "Accessibility%": "{:.1f}"}),
            use_container_width=True
        )

        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)
        st.download_button("💾 Download Full Analysis", data=csv_buffer.getvalue(), file_name=f"{seq_name}_Analysis.csv", use_container_width=True)
