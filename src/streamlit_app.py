import streamlit as st

def check_password():
    """Simple password check."""
    password = st.sidebar.text_input("Enter a password:", type="password")
    if password == "SnowMan0122":
        return True
    else:
        return False

if check_password():
    st.title("My Secure App")
    st.write("Here's our first attempt at using data to create a table:")
    st.write(pd.DataFrame({
        'first column': [1, 2, 3, 4],
        'second column': [10, 20, 30, 40]
    }))
else:
    st.sidebar.error("The password you entered is incorrect.")


