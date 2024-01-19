import streamlit as st
import io
import os
import tempfile
from func1_blastORF import func1_blastORF

def check_password():
	password = st.sidebar.text_input("Enter a password:", type="password")
	if password == "SnowMan0122":
		return True
	else:
		return False

def blastp():
	uploaded_file = st.file_uploader("Upload file", type=['fasta','fa'])
	output_path = tempfile.mkdtemp() 

	if uploaded_file is not None:
		input_file_path = os.path.join(output_path, uploaded_file.name)
		with open(input_file_path, "wb") as f:
			f.write(uploaded_file.getbuffer())
		
		error = func1_blastORF(input_file_path, output_path)
		if error:
			st.error(f"BLASTp error: {error}")
		else:
			st.success("BLASTp analysis completed.")

		if st.button("Download BLASTp Result"):
			output_file = os.path.join(output_path, "blastp_output.txt")
			with open(output_file, "rb") as f:
				data = f.read()
				st.download_button(label="Download BLASTp Output", data=data, file_name="blastp_output.txt", mime="text/plain")
			


if check_password():
	st.title("HERVOminer")
	st.write("test")
	blastp()
	

else:
	st.sidebar.error("The password you entered is incorrect.")

