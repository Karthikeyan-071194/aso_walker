import streamlit as st
import pandas as pd
import io
import requests
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

def get_vienna_fold_api(sequence):
    api_url = f"https://api.vienna-rna.org/rnafold?seq={sequence}"
    try:
        response = requests.get(api_url, timeout=10)
        if response.status_code == 200:
            return response.json().get('structure')
    except:
        return None

def get_internal_fold(sequence):
    n = len(sequence)
    energies = {'CG': 3.4, 'AU': 0.9, 'GU': 0.1}
    dp = np.zeros((n, n))
    for k in range(4, n):
        for i in range(n - k):
            j = i + k
            res = [dp[i+1][j], dp[i][j-1]]
            pair = "".join(sorted([sequence[i], sequence[j]]))
            if pair in energies: res.append(dp[i+1][j-1] + energies[pair])
            for t in range(i + 1, j): res.append(dp[i][t] + dp[t+1][j])
            dp[i][j] = max(res)
    structure = ["."] * n
    stack = [(0, n - 1)]
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
                    if dp[i][j] == dp[i][k] + dp[k+1][j]:
                        stack.append((k+1, j)); stack.append((i, k))
                        break
    return "".join(structure)

def find_best_match(target, database_seq):
    """Finds best match in a sequence and returns the number of mismatches."""
    t_len = len(target)
    d_len = len(database_seq)
    min_mismatches = t_len
    
    for i in range(d_len - t_len + 1):
        window = database_seq[i : i + t_len]
        mismatches = sum(1 for a, b in zip(target, window) if a != b)
        if mismatches < min_mismatches:
            min_mismatches = mismatches
        if min_mismatches == 0: break
    return min_mismatches

# --- UI SECTION ---

st.markdown("<h1 style='text-align: center;'>🧬 ASO Walker Pro</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; color: gray;'>Multi-Sequence Analysis with Mismatch Tolerance</p>", unsafe_allow_html=True)

# 1. Primary Target Input
st.subheader("1. Primary Target Sequence")
c_m1, c_m2 = st.columns([1, 3])
with c_m1:
    main_name = st.text_input("Main Sequence Title", value="WildType_Target")
with c_m2:
    main_raw = st.text_area("Primary Sequence (Walk is performed here)", placeholder="ATCG...", height=100)
main_clean = "".join(main_raw.upper().split())

# 2. Conservation Database
st.subheader("2. Variant Database")
if 'extra_seqs' not in st.session_state:
    st.session_state.extra_seqs = [{"title": "Variant_1", "seq": ""}]

def add_box(): st.session_state.extra_seqs.append({"title": f"Variant_{len(st.session_state.extra_seqs)+1}", "seq": ""})
def remove_box(): 
    if len(st.session_state.extra_seqs) > 1: st.session_state.extra_seqs.pop()

for i, box in enumerate(st.session_state.extra_seqs):
    c1, c2 = st.columns([1, 3])
    with c1: st.session_state.extra_seqs[i]["title"] = st.text_input(f"Title {i+1}", value=box["title"], key=f"t_{i}")
    with c2: st.session_state.extra_seqs[i]["seq"] = st.text_area(f"Seq {i+1}", value=box["seq"], key=f"s_{i}", height=68)

st.button("➕ Add Variant", on_click=add_box)
if len(st.session_state.extra_seqs) > 1: st.button("➖ Remove Variant", on_click=remove_box)

# 3. Settings
st.divider()
c_s1, c_s2, c_s3, c_s4 = st.columns(4)
with c_s1: aso_size = st.number_input("ASO Size (bp)", min_value=1, value=20)
with c_s2: step_size = st.slider("Step Size", 1, 10, 1)
with c_s3: mm_limit = st.slider("Mismatch Tolerance", 0, 3, 0)
with c_s4: st.caption("0 = Perfect match required. 1-3 = Allows near-matches in variants.")

# 4. Execution
if st.button("Run Advanced Conservation Analysis", type="primary", use_container_width=True):
    if not main_clean:
        st.error("Missing Primary Target Sequence.")
    else:
        with st.spinner("Analyzing structures and variant matches..."):
            dot_bracket = get_vienna_fold_api(main_clean) or get_internal_fold(main_clean)
            db = [{"title": i["title"], "seq": "".join(i["seq"].upper().split())} for i in st.session_state.extra_seqs if i["seq"]]

            results = []
            for i in range(0, len(main_clean) - aso_size + 1, step_size):
                target_site = main_clean[i : i + aso_size]
                aso_seq = get_reverse_complement(target_site)
                window_struct = dot_bracket[i : i + aso_size]
                gc, tm = calculate_metrics(aso_seq)
                
                # Check variants
                hits, fails = [], []
                for entry in db:
                    mismatches = find_best_match(target_site, entry["seq"])
                    if mismatches <= mm_limit:
                        hits.append(f"{entry['title']}({mismatches}mm)")
                    else:
                        fails.append(entry["title"])
                
                status = "CONSERVED" if not fails else f"VARIABLE ({len(hits)}/{len(db)})"
                
                results.append({
                    "ASO_ID": f"{main_name}_{len(results) + 1}",
                    "Region": f"{i+1}-{i+aso_size}",
                    "ASO_Sequence": aso_seq,
                    "GC%": gc, "Tm_C": tm,
                    "Accessibility%": round((window_struct.count(".") / aso_size) * 100, 1),
                    "Status": status,
                    "Matched_In": ", ".join(hits),
                    "Missing_In": ", ".join(fails) if fails else "None"
                })

            df = pd.DataFrame(results)

            # High-Contrast Alignment Map
            st.subheader("Target Structure Alignment")
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

            st.subheader("Analysis Results")
            st.dataframe(df.style.background_gradient(subset=['Accessibility%'], cmap='Greens'), use_container_width=True)

            csv_buf = io.StringIO()
            df.to_csv(csv_buf, index=False)
            st.download_button("💾 Download Results", data=csv_buf.getvalue(), file_name=f"{main_name}_Analysis.csv", use_container_width=True)

with st.sidebar:
    st.info(f"**Tolerance:** Allows up to {mm_limit} mismatches in variants.")
    st.divider()
    st.markdown("**ASO Walker Pro**")
