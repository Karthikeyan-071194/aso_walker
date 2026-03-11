import streamlit as st
import pandas as pd
import io
import requests
import numpy as np

# 1. Page Configuration
st.set_page_config(page_title="ASO Walker Pro", page_icon="🧬", layout="wide")

# --- SCORING MATRICES DATA ---
# Values extracted from user-provided images image_eef1da.png
MOE_MATRIX = { # 5-10-5 MOE, all PS
    'A': {1: -3.2, 2: -3.9, 3: -3.6, 4: -2.3, 5: -1.4, 6: -1.4, 7: -3.1, 8: -2.9, 9: -2.2, 10: -2.0, 13: -2.3, 14: -3.0, 15: -2.1, 16: -2.6, 17: -3.8, 18: -4.2, 19: -4.5, 20: -4.1},
    'C': {1: 1.9, 2: 3.1, 3: 4.0, 4: 2.5, 5: 1.9, 13: 2.1, 14: 1.7, 15: 4.2, 16: 4.6, 17: 3.0, 18: 3.8, 19: 4.7},
    'G': {1: 3.1, 2: 1.5, 13: 1.5, 15: -1.9, 17: 1.4, 18: 1.8, 20: -1.6},
    'T': {1: -1.9, 7: 2.4, 8: 2.4, 9: 2.8, 10: 1.9, 11: 1.5}
}

CET_MATRIX = { # 3-10-3 cEt, all PS
    'A': {2: -1.1, 3: 1.1, 4: 2.3, 5: -4.9, 9: 1.1, 13: -2.4, 14: -3.2, 15: -3.1},
    'C': {1: -3.0, 2: -2.6, 3: -4.1, 4: -7.1, 5: -4.0, 6: -4.7, 7: -3.6, 8: -4.5, 9: -4.7, 10: -2.4, 11: -2.4, 12: -1.5, 13: 1.2, 14: 1.1},
    'G': {7: -2.9, 8: -2.0, 11: 1.3, 12: 1.1, 14: -1.9, 15: 2.4},
    'T': {1: 1.7, 2: 4.5, 3: 3.4, 4: 4.3, 5: 8.4, 6: 7.0, 7: 5.0, 8: 4.5, 9: 2.7, 12: 2.0, 13: 2.9, 14: 1.7}
}

def calculate_mod_score(sequence, matrix):
    score = 0.0
    for i, base in enumerate(sequence):
        pos = i + 1
        if base in matrix and pos in matrix[base]:
            score += matrix[base][pos]
    return round(score, 2)

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

def get_vienna_fold_api(sequence):
    api_url = f"https://api.vienna-rna.org/rnafold?seq={sequence}"
    try:
        response = requests.get(api_url, timeout=10)
        if response.status_code == 200:
            return response.json().get('structure')
    except:
        return None

# --- UI SECTION ---

st.markdown("<h1 style='text-align: center;'>🧬 ASO Walker Pro</h1>", unsafe_allow_html=True)

with st.sidebar:
    st.header("Chemistry Settings")
    mod_choice = st.selectbox(
        "Select ASO Modification",
        ["5-10-5 MOE, all PS", "3-10-3 cEt, all PS"]
    )
    current_matrix = MOE_MATRIX if "MOE" in mod_choice else CET_MATRIX
    st.info(f"Scoring applied based on the {mod_choice} efficacy matrix.")
    st.divider()
    mm_limit = st.slider("Mismatch Tolerance", 0, 3, 0)
    st.markdown("**ASO Walker Pro: Predictive Efficacy**")

# Inputs
st.subheader("1. Primary Target Sequence")
c_m1, c_m2 = st.columns([1, 3])
with c_m1:
    main_name = st.text_input("Project Name", value="WT_Target")
with c_m2:
    main_raw = st.text_area("Primary Sequence", placeholder="ATCG...", height=100)
main_clean = "".join(main_raw.upper().split())

# Dynamic Variants
st.subheader("2. Variant Database")
if 'extra_seqs' not in st.session_state:
    st.session_state.extra_seqs = [{"title": "Variant_1", "seq": ""}]

def add_box(): st.session_state.extra_seqs.append({"title": f"Variant_{len(st.session_state.extra_seqs)+1}", "seq": ""})
for i, box in enumerate(st.session_state.extra_seqs):
    c1, c2 = st.columns([1, 3])
    with c1: st.session_state.extra_seqs[i]["title"] = st.text_input(f"Title {i+1}", value=box["title"], key=f"t_{i}")
    with c2: st.session_state.extra_seqs[i]["seq"] = st.text_area(f"Seq {i+1}", value=box["seq"], key=f"s_{i}", height=68)
st.button("➕ Add Variant", on_click=add_box)

st.divider()
c_s1, c_s2 = st.columns(2)
with c_s1: aso_size = st.number_input("ASO Size (bp)", min_value=1, value=20 if "MOE" in mod_choice else 16)
with c_s2: step_size = st.slider("Step Size", 1, 10, 1)

if st.button("Generate Predictive Analysis", type="primary", use_container_width=True):
    if not main_clean:
        st.error("Missing Target Sequence.")
    else:
        with st.spinner("Analyzing structure and predicting efficacy..."):
            dot_bracket = get_vienna_fold_api(main_clean) or ("." * len(main_clean))
            db = [{"title": i["title"], "seq": "".join(i["seq"].upper().split())} for i in st.session_state.extra_seqs if i["seq"]]

            results = []
            for i in range(0, len(main_clean) - aso_size + 1, step_size):
                target_site = main_clean[i : i + aso_size]
                aso_seq = get_reverse_complement(target_site)
                window_struct = dot_bracket[i : i + aso_size]
                gc, tm = calculate_metrics(aso_seq)
                
                # Modification Score from Matrix
                mod_score = calculate_mod_score(aso_seq, current_matrix)
                
                results.append({
                    "ASO_ID": f"{main_name}_{len(results) + 1}",
                    "Region": f"bp_{i+1}_to_{i+aso_size}",
                    "ASO_Sequence": aso_seq,
                    "GC%": gc, "Tm_C": tm,
                    "Accessibility%": round((window_struct.count(".") / aso_size) * 100, 1),
                    "Mod_Efficacy_Score": mod_score,
                    "Chemistry": mod_choice
                })

            df = pd.DataFrame(results)

            # Alignment Map
            r_list = [" "] * len(main_clean)
            for j in range(len(main_clean)):
                if (j+1) == 1 or (j+1) % 10 == 0:
                    for k, d in enumerate(str(j+1)):
                        if j+k < len(main_clean): r_list[j+k] = d
            
            st.markdown(f"""
                <div style="overflow-x: auto; white-space: pre; font-family: 'Courier New', monospace; 
                            background-color: #1E1E1E; padding: 25px; border-radius: 10px; line-height: 2.0; color: white;">
<span style="color: #FF5F5F;">NUM:</span> {"".join(r_list)}
<span style="color: #569CD6;">SEQ:</span> {main_clean}
<span style="color: #4EC9B0;">STR:</span> {dot_bracket}
                </div>""", unsafe_allow_html=True)

            st.subheader("Predictive Results")
            # Color code score: Blue for positive (better), Red for negative
            st.dataframe(df.style.background_gradient(subset=['Mod_Efficacy_Score'], cmap='RdYlBu'), use_container_width=True)

            csv_buf = io.StringIO()
            df.to_csv(csv_buf, index=False)
            st.download_button("💾 Download Predictive CSV", data=csv_buf.getvalue(), file_name=f"{main_name}_Predictive_Analysis.csv", use_container_width=True)




