import streamlit as st

st.set_page_config(
    page_title="Trial Data Dashboard",
    page_icon="🧪",
    layout="wide",
)

st.title("Trial Data Dashboard")
st.markdown(
    """
    This dashboard presents results from the clinical trial analysis pipeline.

    Use the tabs below to explore:
    - the overall project summary
    - Part 2 frequency results
    - Part 3 statistical analysis
    - Part 4 baseline subset analysis
    """
)

tab_overview, tab_part2, tab_part3, tab_part4 = st.tabs(
    [
        "Overview",
        "Part 2: Frequency Summary",
        "Part 3: Statistical Analysis",
        "Part 4: Baseline Subset",
    ]
)

with tab_overview:
    st.header("Overview")
    st.info("Dashboard skeleton is working.")
    st.write("This tab will later summarize the generated outputs and project structure.")

with tab_part2:
    st.header("Part 2: Frequency Summary")
    st.info("Part 2 content will be added next.")
    st.write("This section will display the sample-level population frequency summary table.")

with tab_part3:
    st.header("Part 3: Statistical Analysis")
    st.info("Part 3 content will be added next.")
    st.write("This section will display GEE results and responder vs non-responder plots.")

with tab_part4:
    st.header("Part 4: Baseline Subset")
    st.info("Part 4 content will be added next.")
    st.write("This section will display the baseline melanoma PBMC miraclib subset summaries.")
