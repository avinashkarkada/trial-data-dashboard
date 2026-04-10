from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="Trial Data Dashboard",
    page_icon="🧪",
    layout="wide",
)

OUTPUT_DIR = Path("outputs")

PART2_FILE = OUTPUT_DIR / "sample_population_frequencies.csv"
PART3_FILTERED_FILE = OUTPUT_DIR / "part3_filtered_frequencies.csv"
PART3_STATS_FILE = OUTPUT_DIR / "part3_population_statistics_gee.csv"
PART3_SUMMARY_FILE = OUTPUT_DIR / "part3_gee_model_summaries.txt"
PART3_PLOT_FILE = OUTPUT_DIR / "part3_boxplots.png"

PART4_BASELINE_FILE = OUTPUT_DIR / "part4_baseline_samples.csv"
PART4_PROJECT_FILE = OUTPUT_DIR / "part4_samples_per_project.csv"
PART4_RESPONSE_FILE = OUTPUT_DIR / "part4_subjects_by_response.csv"
PART4_SEX_FILE = OUTPUT_DIR / "part4_subjects_by_sex.csv"
PART4_SUMMARY_FILE = OUTPUT_DIR / "part4_summary.txt"


def load_csv_if_exists(path: Path):
    if path.exists():
        return pd.read_csv(path)
    return None


def load_text_if_exists(path: Path):
    if path.exists():
        return path.read_text(encoding="utf-8")
    return None


part2_df = load_csv_if_exists(PART2_FILE)
part3_filtered_df = load_csv_if_exists(PART3_FILTERED_FILE)
part3_stats_df = load_csv_if_exists(PART3_STATS_FILE)
part3_summary_text = load_text_if_exists(PART3_SUMMARY_FILE)
part4_baseline_df = load_csv_if_exists(PART4_BASELINE_FILE)
part4_project_df = load_csv_if_exists(PART4_PROJECT_FILE)
part4_response_df = load_csv_if_exists(PART4_RESPONSE_FILE)
part4_sex_df = load_csv_if_exists(PART4_SEX_FILE)
part4_summary_text = load_text_if_exists(PART4_SUMMARY_FILE)

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

    if not OUTPUT_DIR.exists():
        st.error("The outputs directory was not found. Run the analysis pipeline first.")
    else:
        st.success(f"Found outputs directory: {OUTPUT_DIR}")

    file_status = {
        "Part 2 frequency summary": PART2_FILE.exists(),
        "Part 3 filtered frequencies": PART3_FILTERED_FILE.exists(),
        "Part 3 GEE statistics": PART3_STATS_FILE.exists(),
        "Part 3 GEE summaries": PART3_SUMMARY_FILE.exists(),
        "Part 3 boxplot": PART3_PLOT_FILE.exists(),
        "Part 4 baseline samples": PART4_BASELINE_FILE.exists(),
        "Part 4 samples per project": PART4_PROJECT_FILE.exists(),
        "Part 4 subjects by response": PART4_RESPONSE_FILE.exists(),
        "Part 4 subjects by sex": PART4_SEX_FILE.exists(),
        "Part 4 summary": PART4_SUMMARY_FILE.exists(),
    }

    status_df = pd.DataFrame(
        {
            "Output": list(file_status.keys()),
            "Found": ["Yes" if found else "No" for found in file_status.values()],
        }
    )
    st.dataframe(status_df, use_container_width=True, hide_index=True)

    st.subheader("Loaded data preview")
    st.write(
        {
            "part2_rows": 0 if part2_df is None else len(part2_df),
            "part3_filtered_rows": 0 if part3_filtered_df is None else len(part3_filtered_df),
            "part3_stats_rows": 0 if part3_stats_df is None else len(part3_stats_df),
            "part4_baseline_rows": 0 if part4_baseline_df is None else len(part4_baseline_df),
        }
    )

with tab_part2:
    st.header("Part 2: Frequency Summary")

    if part2_df is None:
        st.warning("Part 2 output file not found. Run the pipeline first.")
    else:
        st.markdown(
            """
            This table shows the relative frequency of each immune cell population in each sample.
            For every sample, the total cell count is computed across the five populations, and each
            population's frequency is reported as a percentage of that total.
            """
        )

        col1, col2, col3 = st.columns(3)
        col1.metric("Rows", f"{len(part2_df):,}")
        col2.metric("Unique samples", f"{part2_df['sample'].nunique():,}")
        col3.metric("Populations", f"{part2_df['population'].nunique()}")

        st.subheader("Filters")

        available_populations = sorted(part2_df["population"].dropna().unique().tolist())
        selected_populations = st.multiselect(
            "Select population(s)",
            options=available_populations,
            default=available_populations,
        )

        sample_search = st.text_input(
            "Search sample ID",
            placeholder="Type part of a sample ID...",
        )

        max_rows = st.slider(
            "Maximum rows to display",
            min_value=10,
            max_value=500,
            value=100,
            step=10,
        )

        filtered_part2_df = part2_df.copy()

        if selected_populations:
            filtered_part2_df = filtered_part2_df[
                filtered_part2_df["population"].isin(selected_populations)
            ]

        if sample_search.strip():
            filtered_part2_df = filtered_part2_df[
                filtered_part2_df["sample"].str.contains(sample_search, case=False, na=False)
            ]

        st.subheader("Filtered frequency table")
        st.write(f"Showing {min(len(filtered_part2_df), max_rows):,} of {len(filtered_part2_df):,} rows")

        st.dataframe(
            filtered_part2_df.head(max_rows),
            use_container_width=True,
            hide_index=True,
        )

        with st.expander("Column descriptions"):
            st.markdown(
                """
                - **sample**: sample identifier
                - **total_count**: total number of cells across all five populations in that sample
                - **population**: immune cell population name
                - **count**: raw cell count for that population in the sample
                - **percentage**: relative frequency of that population within the sample
                """
            )

with tab_part3:
    st.header("Part 3: Statistical Analysis")
    if part3_stats_df is None:
        st.warning("Part 3 GEE output file not found. Run the pipeline first.")
    else:
        st.success("Part 3 output files loaded successfully.")
        st.write("Detailed display will be added next.")

with tab_part4:
    st.header("Part 4: Baseline Subset")
    if part4_baseline_df is None:
        st.warning("Part 4 output files not found. Run the pipeline first.")
    else:
        st.success("Part 4 output files loaded successfully.")
        st.write("Detailed display will be added next.")
