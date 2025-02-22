import streamlit as st
import pandas as pd
import plotly.express as px
import io

def load_example_data():
    """Load the example dataset and store it in session state."""
    st.session_state.df = px.data.tips()
    st.session_state.data_loaded = True
    st.success(f"Loaded example dataset with {len(st.session_state.df)} rows")

def detect_delimiter(content):
    """Detects the delimiter in text-based files."""
    if content.count(',') > content.count('\t') and content.count(',') > content.count(' '):
        return ','
    elif content.count('\t') > content.count(' '):
        return '\t'
    else:
        return ' '

def load_uploaded_file(uploaded_file):
    """Loads uploaded file into a Pandas DataFrame."""
    try:
        file_ext = uploaded_file.name.split('.')[-1].lower()
        
        if file_ext in ["csv", "tsv", "txt"]:
            content = uploaded_file.getvalue().decode()
            delimiter = detect_delimiter(content)
            df = pd.read_csv(uploaded_file, delimiter=delimiter)
        elif file_ext == "xlsx":
            df = pd.read_excel(uploaded_file)
        elif file_ext == "json":
            df = pd.read_json(uploaded_file)
        elif file_ext == "parquet":
            df = pd.read_parquet(uploaded_file)
        else:
            st.error("Unsupported file format. Please upload CSV, TSV, TXT, Excel, JSON, or Parquet.")
            return
        
        st.session_state.df = df
        st.session_state.data_loaded = True
        st.success(f"Successfully loaded {uploaded_file.name} with {len(df)} rows")
    except Exception as e:
        st.error(f"Could not read file: {str(e)}")
        st.session_state.data_loaded = False

def load_manual_data(manual_text):
    """Processes manually entered tabular data into a DataFrame."""
    try:
        delimiter = detect_delimiter(manual_text)
        df = pd.read_csv(io.StringIO(manual_text), delimiter=delimiter)
        st.session_state.df = df
        st.session_state.data_loaded = True
        st.success("Successfully loaded manual data.")
    except Exception as e:
        st.error(f"Could not process manual data: {str(e)}")
        st.session_state.data_loaded = False

def main():
    st.title("Interactive Data Visualization App")

    # Initialize session state variables
    if "df" not in st.session_state:
        st.session_state.df = None
        st.session_state.data_loaded = False

    # File Upload Section
    uploaded_file = st.file_uploader(
        "Upload a file (CSV, TSV, TXT, Excel, JSON, Parquet)",
        type=["csv", "tsv", "txt", "xlsx", "json", "parquet"]
    )
    if uploaded_file is not None:
        load_uploaded_file(uploaded_file)

    # Manual Data Input Section
    manual_text = st.text_area(
        "Or enter/paste your own data (CSV, TSV, or space-separated values)",
        height=150,
        help="Ensure headers are included in the first row"
    )
    if st.button("Load manual data"):
        if manual_text.strip():
            load_manual_data(manual_text)
        else:
            st.warning("Please enter some data before clicking the button.")

    # Button to load example data
    if st.button("Use example data (Tips dataset)"):
        load_example_data()

    # Prevent proceeding if no data is loaded
    if not st.session_state.data_loaded:
        st.info("Please upload a file, enter data manually, or use the example dataset.")
        return

    df = st.session_state.df

    if st.checkbox("Show raw data"):
        st.dataframe(df.head(10))

    st.sidebar.header("Plot Configuration")
    plot_type = st.sidebar.selectbox(
        "Choose plot type",
        ["Scatter", "Line", "Bar", "Histogram", "Box", "Pair Plot"]
    )

    numeric_cols = df.select_dtypes(include='number').columns.tolist()
    categorical_cols = df.select_dtypes(exclude='number').columns.tolist()
    all_cols = df.columns.tolist()

    if plot_type in ["Scatter", "Line", "Bar"]:
        x_axis = st.sidebar.selectbox("X Axis", all_cols, index=0)
        y_axis = st.sidebar.selectbox("Y Axis", numeric_cols, index=0)
    elif plot_type in ["Box", "Histogram"]:
        x_axis = st.sidebar.selectbox("X Axis", categorical_cols if plot_type == "Box" else all_cols)
        y_axis = st.sidebar.selectbox("Y Axis", numeric_cols) if plot_type == "Box" else None
    else:  # Pair plot
        dimensions = st.sidebar.multiselect("Select dimensions", numeric_cols, default=numeric_cols[:3])

    # Color selection for applicable plot types
    color_var = None
    if plot_type in ["Scatter", "Bar"]:
        color_var = st.sidebar.selectbox("Color by", ["None"] + all_cols, index=0)
        color_var = None if color_var == "None" else color_var

    st.subheader(f"{plot_type} Plot")
    try:
        if plot_type == "Scatter":
            fig = px.scatter(df, x=x_axis, y=y_axis, color=color_var)
        elif plot_type == "Line":
            fig = px.line(df, x=x_axis, y=y_axis)
        elif plot_type == "Bar":
            fig = px.bar(df, x=x_axis, y=y_axis, color=color_var)
        elif plot_type == "Histogram":
            fig = px.histogram(df, x=x_axis)
        elif plot_type == "Box":
            fig = px.box(df, x=x_axis, y=y_axis)
        elif plot_type == "Pair Plot":
            fig = px.scatter_matrix(df[dimensions])

        st.plotly_chart(fig, use_container_width=True)
    except Exception as e:
        st.error(f"Error generating plot: {str(e)}")

    # Footer
    st.markdown("---")
    st.markdown("2025 | Michael Bohl")

if __name__ == "__main__":
    main()
