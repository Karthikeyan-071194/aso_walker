import streamlit as st
import pandas as pd
import io
import requests
import numpy as np

# 1. Page Configuration
st.set_page_config(page_title="ASO Walker Pro", page_icon="🧬", layout="wide")

# --- MOTIF DATA ---
LOW_EFFICACY_MOTIFS = {'GGGG': '24.0%', 'AAAA': '29.0%', 'TAAA': '32.0%', 'CTAA': '34.0%', 'CCTA': '35.0%'}
HIGH_EFFICACY_MOTIFS = {'TTGT': '53.0%', 'GTAT': '54.0%', 'CGTA': '54.0%', 'GTCG': '54.7%', 'GCGT': '57.0%'}

# --- SCORING MATRICES DATA ---
MOE_MATRIX = { 
    'A': {1: -3.2, 2: -3.9, 3: -3.6, 4: -2.3, 5: -1.4, 6: -1.4, 7: -3.1, 8: -2.9, 9: -2.2, 10: -2.0, 13: -2.3, 14: -3.0, 15: -2.1, 16: -2.6, 17: -3.8, 18: -4.2, 19: -4.5, 20: -4.1},
    'C': {1: 1.9, 2: 3.1, 3: 4.0, 4: 2.5, 5: 1.9, 13: 2.1, 14: 1.7, 15: 4.2, 16: 4.6, 17: 3.0, 18: 3.8, 19: 4.7},
    'G': {1: 3.1, 2: 1.5, 13: 1.5, 15: -1.9, 17: 1.4, 18: 1.8, 20: -1.6},
    'T': {1: -1.9, 7: 2.4, 8: 2.4, 9: 2.8, 10: 1.9, 11: 1.5}
}

CET_MATRIX = { 
    'A': {2: -1.1, 3: 1.1, 4: 2.3, 5: -4.9, 9: 1.1, 13: -2.4, 14: -3.2, 15: -3.1},
    'C': {1: -3.0, 2: -2.6, 3: -4.1, 4: -7.1, 5: -4.0, 6: -4.7, 7: -3.6, 8: -4.5, 9: -4.7, 10: -2.4, 11: -2.4, 12: -1.5, 13: 1.2, 14: 1.1},
    'G': {7: -2.9, 8: -2.0, 11: 1.3, 12: 1.1, 14: -1.9, 15: 2.4},
    'T': {1: 1.7, 2: 4.5, 3: 3.4, 4: 4.3, 5: 8.4, 6: 7.0, 7: 5.0, 8: 4.5, 9: 2.7, 12: 2.0, 13: 2.9, 14: 1.7}
}

def calculate_mod_score(sequence, matrix):
    score = sum(matrix.get(base, {}).get(i+1, 0.0) for i, base in enumerate(sequence))
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
    return round(gc_cont, 2), tm

def get_vienna_fold_api(sequence):
    api_url = f"https://api.vienna-rna.org/rnafold?seq={sequence}"
    try:
        response = requests.get(api_url, timeout=10)
        return response.json().get('structure') if response.status_code == 200 else None
    except: return None

def get_internal_fold(sequence):
    n = len(sequence)
    energies, dp = {'CG': 3.4, 'AU': 0.9, 'GU': 0.1}, np.zeros((n, n))
    for k in range(4, n):
        for i in range(n - k):
            j = i + k
            res = [dp[i+1][j], dp[i][j-1]]
            pair = "".join(sorted([sequence[i], sequence[j]]))
            if pair in energies: res.append(dp[i+1][j-1] + energies[pair])
            for t in range(i + 1, j): res.append(dp[i][t] + dp[t+1][j])
            dp[i][j] = max(res)
    structure, stack = ["."] * n, [(0, n - 1)]
    while stack:
        i, j = stack.pop()
        if i >= j - 3: continue 
        elif dp[i][j] == dp[i+1][j]: stack.append((i+1, j))
        elif dp[i][j] == dp[i][j-1]: stack.append((i, j-1))
        else:
            pair = "".join(sorted([sequence[i], sequence[j]]))
            if pair in energies and dp[i][j] == dp[i+1][j-1] + energies[pair]:
                structure[i], structure[j] = "(", ")"
                stack.append((i+1, j-1))
            else:
                for k in range(i+1, j):
                    if dp[i][j] == dp[i][k] + dp[t+1][j]:
                        stack.append((k+1, j)); stack.append((i, k))
                        break
    return "".join(structure)

def find_best_match(target, database_seq):
    t_len, d_len, min_mm = len(target), len(database_seq), len(target)
    for i in range(d_len - t_len + 1):
        mm = sum(1 for a, b in zip(target, database_seq[i:i+t_len]) if a != b)
        min_mm = min(min_mm, mm)
        if min_mm == 0: break
    return min_mm

# --- UI SECTION ---
st.markdown("<h1 style='text-align: center;'>🧬 ASO Walker Pro</h1>", unsafe_allow_html=True)

with st.sidebar:
    st.header("Settings")
    mod_choice = st.selectbox("Chemistry Choice", ["5-10-5 MOE, all PS", "3-10-3 cEt, all PS"])
    mm_limit = st.slider("Mismatch Tolerance", 0, 3, 0)
    current_matrix = MOE_MATRIX if "MOE" in mod_choice else CET_MATRIX

raw_seq = st.text_area("Primary Target Sequence", placeholder="Paste sequence here...", height=100)
clean_seq = "".join(raw_seq.upper().split())

col1, col2, col3 = st.columns(3)
with col1: seq_name = st.text_input("Project Title", value="Target_1")
with col2: aso_size = st.number_input("ASO Size (bp)", min_value=1, value=20 if "MOE" in mod_choice else 16)
with col3: step_size = st.slider("Stride (Step)", 1, 10, 1)

st.subheader("Variant Database")
if 'extra_seqs' not in st.session_state: st.session_state.extra_seqs = [{"title": "Variant_1", "seq": ""}]
def add_box(): st.session_state.extra_seqs.append({"title": f"Variant_{len(st.session_state.extra_seqs)+1}", "seq": ""})
def remove_box(): 
    if len(st.session_state.extra_seqs) > 1: st.session_state.extra_seqs.pop()

for i, box in enumerate(st.session_state.extra_seqs):
    c1, c2 = st.columns([1, 3])
    with c1: st.session_state.extra_seqs[i]["title"] = st.text_input(f"Title {i+1}", value=box["title"], key=f"t_{i}")
    with c2: st.session_state.extra_seqs[i]["seq"] = st.text_area(f"Seq {i+1}", value=box["seq"], key=f"s_{i}", height=68)

st.button("➕ Add Variant", on_click=add_box)
if len(st.session_state.extra_seqs) > 1: st.button("➖ Remove Variant", on_click=remove_box)

if st.button("Generate Complete Analysis", type="primary", use_container_width=True):
    if not clean_seq: st.error("Please enter a target sequence.")
    else:
        with st.spinner("Analyzing variants and folding..."):
            dot_bracket = get_vienna_fold_api(clean_seq) or get_internal_fold(clean_seq)
            db = [{"title": i["title"], "seq": "".join(i["seq"].upper().split())} for i in st.session_state.extra_seqs if i["seq"]]

        r_list = [" "] * len(clean_seq)
        for i in range(len(clean_seq)):
            if (i+1)==1 or (i+1)%10==0:
                for j, d in enumerate(str(i+1)):
                    if i+j < len(clean_seq): r_list[i+j] = d
        
        st.subheader("Thermodynamic Alignment Map")
        st.markdown(f"""<div style="overflow-x: auto; white-space: pre; font-family: 'Courier New', monospace; background-color: #1E1E1E; padding: 25px; border-radius: 10px; line-height: 2.0; border: 1px solid #333;">
<span style="color: #FF5F5F; font-weight: bold;">NUM:</span> <span style="color: #FFFFFF;">{"".join(r_list)}</span>
<span style="color: #569CD6; font-weight: bold;">SEQ:</span> <span style="color: #DCDCAA; font-weight: bold;">{clean_seq}</span>
<span style="color: #4EC9B0; font-weight: bold;">STR:</span> <span style="color: #CE9178;">{dot_bracket}</span></div>""", unsafe_allow_html=True)

        results = []
        for i in range(0, len(clean_seq) - aso_size + 1, step_size):
            target_site = clean_seq[i : i + aso_size]
            aso_seq = get_reverse_complement(target_site)
            window_struct = dot_bracket[i : i + aso_size]
            gc, tm = calculate_metrics(aso_seq)
            
            # Restoration of Conservation Search
            hits, fails = [], []
            for entry in db:
                mm = find_best_match(target_site, entry["seq"])
                if mm <= mm_limit: hits.append(f"{entry['title']}({mm}mm)")
                else: fails.append(entry["title"])
            
            mod_score = calculate_mod_score(aso_seq, current_matrix)

            # Motif Logic
            found_motifs = []
            font_color = "#DCDCAA" # Default
            for m, val in LOW_EFFICACY_MOTIFS.items():
                if m in aso_seq: 
                    found_motifs.append(f"{m}({val})")
                    font_color = "red"
            for m, val in HIGH_EFFICACY_MOTIFS.items():
                if m in aso_seq: 
                    found_motifs.append(f"{m}({val})")
                    font_color = "green"

            results.append({
                "ASO_ID": f"{seq_name}_{len(results) + 1}",
                "Region": f"Bases {i+1} to {i+aso_size}",
                "ASO_Sequence": aso_seq,
                "GC%": gc, "Tm_C": tm,
                "Accessibility%": round((window_struct.count(".") / aso_size) * 100, 2),
                "Mod_Score": mod_score,
                "Motifs": ", ".join(found_motifs) if found_motifs else "None",
                "Status": "CONSERVED" if not fails else f"VAR ({len(hits)}/{len(db)})",
                "Matched_In": ", ".join(hits),
                "Missing_In": ", ".join(fails) if fails else "None"
            })
        
        df = pd.DataFrame(results)

        def color_aso(row):
            return [f'color: {row["Color"]}' if col == 'ASO_Sequence' else '' for col in df.columns]

        st.dataframe(
            df.style.apply(color_aso, axis=1)
            .background_gradient(subset=['Mod_Score'], cmap='RdYlBu')
            .format({"GC%": "{:.2f}", "Accessibility%": "{:.2f}", "Mod_Score": "{:.2f}"}), 
            use_container_width=True
        )

        csv_buffer = io.StringIO()
        df.drop(columns=['Color']).to_csv(csv_buffer, index=False)
        st.download_button("💾 Download Restored CSV", data=csv_buffer.getvalue(), file_name=f"{seq_name}_Analysis.csv", use_container_width=True)
