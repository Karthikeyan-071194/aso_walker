import streamlit as st
import pandas as pd
import io
import re

# 1. Page Configuration
st.set_page_config(page_title="ASO Walker Pro", page_icon="🧬", layout="centered")

# 2. Advanced Biological Functions
def get_reverse_complement(seq):
    trans = str.maketrans('ATCGUatcgu', 'TAGCAtagca')
    return seq.translate(trans)[::-1]

def calculate_metrics(seq):
    """Calculates GC content and estimated Melting Temperature."""
    seq = seq.upper()
    g = seq.count('G')
    c = seq.count('C')
    a = seq.count('A')
    t = seq.count('T')
    u = seq.count('U')
    
    gc_cont = ((g + c) / len(seq)) * 100
    # Basic Wallace Rule for Tm
    tm = 2 * (a + t + u) + 4 * (g + c)
    return round(gc_cont, 1), tm

# 3. Centered Header
st.markdown(
    """
    <div style="text-align: center; padding-bottom: 20px;">
        <img src="https://github.com/Karthikeyan-071194/aso_walker/blob/main/the_walker.png?raw=true" 
             width="150" style="display: block; margin: auto; margin-bottom: 10px;">
        <h1 style="margin: 0; font-family: 'Courier New';">ASO Walker Pro</h1>
        <p style="color: gray; font-family: 'Courier New';">Advanced Screening & Quality Control</p>
    </div>
    """, unsafe_allow_html=True
)

st.divider()

# 4. Inputs
seq_name = st.text_input("Sequence Name", placeholder="e.g., Exon_4_Target")
raw_seq = st.text_area("Target Sequence (5'-3')", placeholder="Paste sequence here...", height=150)

# Cleaning and Validation
clean_seq = "".join(raw_seq.upper().split())
invalid_chars = re.findall(r'[^ATCGU]', clean_seq)

if invalid_chars:
    st.warning(f"⚠️ Warning: Detected non-nucleotide characters: {set(invalid_chars)}")

col1, col2 = st.columns(2)
with col1:
    aso_size = st.number_input("ASO Size (bp)", min_value=1, value=20)
with col2:
    step_size = st.slider("Step Size", 1, 10, 1)

# 5. Execution
if st.button("Generate & Analyze ASOs", type="primary", use_container_width=True):
    if not clean_seq:
        st.error("Missing target sequence.")
    elif len(clean_seq) < aso_size:
        st.error("ASO size exceeds sequence length.")
    else:
        results = []
        for i in range(0, len(clean_seq) - aso_size + 1, step_size):
            window = clean_seq[i : i + aso_size]
            aso_seq = get_reverse_complement(window)
            gc, tm = calculate_metrics(aso_seq)
            
            results.append({
                "Name": f"{seq_name}_{len(results) + 1}",
                "Sequence 5'-3'": aso_seq,
                "Length": len(aso_seq),
                "GC Content (%)": gc,
                "Est. Tm (°C)": tm
            })
        
        df = pd.DataFrame(results)
        
        # --- FIX: FORCING 1 DECIMAL PLACE IN DISPLAY ---
        # We use .format() to tell Streamlit EXACTLY how to show the numbers
        styled_df = df.style.format({
            "GC Content (%)": "{:.1f}", 
            "Est. Tm (°C)": "{:.0f}"
        })
        
        # Display QC Metrics
        st.success(f"Analysis Complete: {len(df)} ASOs generated.")
        
        # Color-coding the GC Content for quick screening
        st.dataframe(df.style.background_gradient(subset=['GC Content (%)'], cmap='RdYlGn', low=0.4, high=0.6), use_container_width=True)
        
        # Export
        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)
        st.download_button("💾 Download Analysis CSV", data=csv_buffer.getvalue(), file_name=f"{seq_name}_Analysis.csv", use_container_width=True)

with st.sidebar:
    st.header("QC Guidelines")
    st.write("**Ideal GC:** 40% - 60%")
    st.write("**Tm Warning:** Low Tm may indicate weak binding.")
    st.divider()
    st.write("Developed for ASO Research")



