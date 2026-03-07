import streamlit as st
import pandas as pd
import io

# 1. Page Configuration (Sets the browser tab title and favicon)
st.set_page_config(
    page_title="ASO Walker",
    page_icon="🧬",
    layout="centered"
)

# 2. Reverse Complement Logic
def get_reverse_complement(seq):
    """Returns the reverse complement of a DNA/RNA sequence."""
    trans = str.maketrans('ATCGUatcgu', 'TAGCAtagca')
    return seq.translate(trans)[::-1]

# 3. Header Section (Logo + Title)
# We create two columns: one narrow for the logo, one wide for the text
col1, col2 = st.columns([4, 4], vertical_alignment="center")

with col1:
    try:
        # Streamlit handles PNGs with transparency natively
        st.image("the_walker.png", width=80)
    except:
        st.write("🧬") # Fallback emoji if image is missing

with col2:
    st.title("ASO Walker")
    st.caption("Generate Antisense Oligonucleotide sequences via sliding window.")

st.divider()

# 4. Input Section
# Streamlit widgets are responsive by default (they stack on mobile)
seq_name = st.text_input("Sequence Name", placeholder="Enter the name of the sequence", help="Base name for your ASOs (e.g., MyGene)")

raw_seq = st.text_area("Target Sequence (5'-3')", placeholder="Paste your DNA/RNA sequence here...", height=200)

aso_size = st.number_input("ASO Size (bp)", min_value=1, max_value=100, value=20, step=1)

# 5. Processing Logic
if st.button("Generate ASO List", type="primary", use_container_width=True):
    # Sanitize input: remove whitespace and convert to upper
    clean_seq = "".join(raw_seq.upper().split())
    
    # Validation
    if not clean_seq:
        st.error("Please enter a target sequence.")
    elif len(clean_seq) < aso_size:
        st.warning(f"Sequence length ({len(clean_seq)}bp) is shorter than ASO size ({aso_size}bp).")
    else:
        # Generate the windows
        results = []
        num_windows = len(clean_seq) - aso_size + 1
        
        for i in range(num_windows):
            target_window = clean_seq[i : i + aso_size]
            aso_seq = get_reverse_complement(target_window)
            results.append({
                "Name": f"{seq_name}_{i+1}",
                "Sequence 5'-3'": aso_seq,
                "Length": len(aso_seq)
            })
        
        # Create DataFrame
        df = pd.DataFrame(results)
        
        # Display Success and Data Preview
        st.success(f"Generated {num_windows} ASO sequences successfully!")
        
        # Show table (Interactive on desktop, scrollable on mobile)
        st.dataframe(df, use_container_width=True)
        
        # 6. Export to CSV
        # We use a buffer so the file is generated in memory
        csv_buffer = io.StringIO()
        df.to_csv(csv_buffer, index=False)
        csv_data = csv_buffer.getvalue()
        
        st.download_button(
            label="💾 Download CSV File",
            data=csv_data,
            file_name=f"{seq_name}_ASO_Results.csv",
            mime="text/csv",
            use_container_width=True
        )

# 7. Sidebar Info
with st.sidebar:
    st.header("How it works")
    st.info("""
    1. The app extracts a window of your specified size.
    2. It moves **1 base at a time** across your sequence.
    3. It calculates the **Reverse Complement** for each window to create the ASO.
    4. The sequence name will be used for naming the ASOs (Example: if sequence name is ABC, the names of the 
    ASOs will be ABC_1, ABC_2,....)
    """)
    st.markdown("---")

    st.markdown("Developed for ASO Research")






