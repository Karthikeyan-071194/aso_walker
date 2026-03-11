# --- ADD THIS AT THE TOP OF YOUR CODE WITH OTHER IMPORTS ---
import io

# --- REPLACE THE DOWNLOAD BUTTON SECTION WITH THIS ---

if st.button("Generate Comprehensive Analysis", type="primary", use_container_width=True):
    # ... (Keep all your existing calculation logic here) ...
    
    # Create an Excel buffer
    excel_buffer = io.BytesIO()
    
    # Use XlsxWriter as the engine to support formatting
    with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
        # Drop the helper 'Color' column before saving
        export_df = df.drop(columns=['Color'])
        export_df.to_excel(writer, index=False, sheet_name='ASO_Analysis')
        
        workbook  = writer.book
        worksheet = writer.sheets['ASO_Analysis']
        
        # Define the formats
        red_font = workbook.add_format({'font_color': 'red'})
        green_font = workbook.add_format({'font_color': 'green'})
        
        # Find the index of the ASO_Sequence column
        col_idx = export_df.columns.get_loc("ASO_Sequence")
        
        # Apply formatting to the Excel file based on the 'Color' logic
        for row_num, color_val in enumerate(df['Color']):
            if color_val == 'red':
                worksheet.write(row_num + 1, col_idx, df.iloc[row_num]['ASO_Sequence'], red_font)
            elif color_val == 'green':
                worksheet.write(row_num + 1, col_idx, df.iloc[row_num]['ASO_Sequence'], green_font)

    st.download_button(
        label="💾 Download Formatted Excel (.xlsx)",
        data=excel_buffer.getvalue(),
        file_name=f"{seq_name}_Analysis.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        use_container_width=True
    )
