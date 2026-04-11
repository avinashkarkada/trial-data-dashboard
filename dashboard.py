from pathlib import Path

import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

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


def format_part3_stats_table(stats_df: pd.DataFrame) -> pd.DataFrame:
    display_df = stats_df.copy()

    numeric_cols = [
        "response_coef",
        "response_p_value",
        "interaction_time7_coef",
        "interaction_time7_p_value",
        "interaction_time14_coef",
        "interaction_time14_p_value",
        "joint_response_time_wald_stat",
        "joint_response_time_p_value",
        "joint_response_time_fdr_bh",
    ]

    for col in numeric_cols:
        if col in display_df.columns:
            display_df[col] = display_df[col].round(6)

    preferred_columns = [
        "population",
        "n_samples",
        "n_subjects",
        "responder_subjects",
        "non_responder_subjects",
        "response_coef",
        "response_p_value",
        "interaction_time7_coef",
        "interaction_time7_p_value",
        "interaction_time14_coef",
        "interaction_time14_p_value",
        "joint_response_time_wald_stat",
        "joint_response_time_df",
        "joint_response_time_p_value",
        "joint_response_time_fdr_bh",
        "significant_fdr_0_05",
    ]

    existing_columns = [col for col in preferred_columns if col in display_df.columns]
    return display_df[existing_columns]


def build_part3_interpretation(stats_df: pd.DataFrame) -> str:
    if stats_df is None or stats_df.empty:
        return "No Part 3 statistics are available."

    if "significant_fdr_0_05" in stats_df.columns:
        significant_df = stats_df[stats_df["significant_fdr_0_05"] == True]
    else:
        significant_df = pd.DataFrame()

    if not significant_df.empty:
        populations = significant_df["population"].tolist()
        return (
            "After Benjamini-Hochberg correction across the five populations, "
            "the following populations show significant response-related differences: "
            + ", ".join(populations)
            + "."
        )

    if "joint_response_time_p_value" in stats_df.columns:
        best_row = stats_df.sort_values("joint_response_time_p_value").iloc[0]
        best_population = best_row["population"]
        best_p = best_row["joint_response_time_p_value"]
        best_fdr = best_row["joint_response_time_fdr_bh"]
        return (
            "No immune cell population remains significant after Benjamini-Hochberg correction "
            f"at FDR < 0.05. The smallest raw joint response/time p-value was observed for "
            f"{best_population} (p = {best_p:.6g}, FDR = {best_fdr:.6g})."
        )

    return "No statistically significant populations were identified."

def build_part3_interactive_boxplot(filtered_df: pd.DataFrame):
    populations = sorted(filtered_df["population"].dropna().unique().tolist())
    response_map = {"no": "Non-responder", "yes": "Responder"}

    fig = make_subplots(
        rows=1,
        cols=len(populations),
        shared_yaxes=True,
        subplot_titles=populations,
    )

    for col_idx, population in enumerate(populations, start=1):
        pop_df = filtered_df[filtered_df["population"] == population].copy()

        for response_value in ["no", "yes"]:
            resp_df = pop_df[pop_df["response"] == response_value].copy()
            response_label = response_map[response_value]

            if resp_df.empty:
                continue

            fig.add_trace(
                go.Box(
                    y=resp_df["percentage"],
                    x=[response_label] * len(resp_df),
                    name=response_label,
                    legendgroup=response_label,
                    showlegend=False,
                    boxpoints="outliers",
                    width=0.55,
                    customdata=[[population, response_label]] * len(resp_df),
                    hovertemplate=(
                        "Population: %{customdata[0]}<br>"
                        "Response: %{customdata[1]}<br>"
                        "Percentage: %{y:.2f}%<extra></extra>"
                    ),
                ),
                row=1,
                col=col_idx,
            )

            q1 = resp_df["percentage"].quantile(0.25)
            median = resp_df["percentage"].median()
            q3 = resp_df["percentage"].quantile(0.75)
            iqr = q3 - q1
            n_samples = len(resp_df)

            fig.add_trace(
                go.Scatter(
                    x=[response_label],
                    y=[median],
                    mode="markers",
                    name=f"{response_label} summary",
                    legendgroup=response_label,
                    showlegend=False,
                    marker={"size": 8, "symbol": "diamond-open"},
                    customdata=[[population, response_label, n_samples, q1, median, q3, iqr]],
                    hovertemplate=(
                        "Population: %{customdata[0]}<br>"
                        "Response: %{customdata[1]}<br>"
                        "Samples: %{customdata[2]}<br>"
                        "Q1: %{customdata[3]:.2f}%<br>"
                        "Median: %{customdata[4]:.2f}%<br>"
                        "Q3: %{customdata[5]:.2f}%<br>"
                        "IQR: %{customdata[6]:.2f}%<extra></extra>"
                    ),
                ),
                row=1,
                col=col_idx,
            )

    fig.update_yaxes(title_text="Relative frequency (%)", row=1, col=1)

    fig.add_annotation(
        text="Response",
        x=0.5,
        y=-0.16,
        xref="paper",
        yref="paper",
        showarrow=False,
    )

    fig.update_layout(
        title="Melanoma PBMC samples treated with miraclib: responder vs non-responder",
        boxmode="group",
        boxgap=0.25,
        boxgroupgap=0.1,
        showlegend=False,
        height=550,
        margin={"l": 40, "r": 20, "t": 80, "b": 80},
    )

    return fig

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

    st.markdown(
        """
        This dashboard was built to help **Bob Loblaw** review immune cell population
        data from a longitudinal clinical trial and share the findings clearly with
        **Yah D'yada**.

        The pipeline begins with the input file cell-counts.csv, loading the data into a 
        relational SQLite database, computing frequencies of immune cells per sample, comparing 
        responder vs. non-responder within a melanoma PBMC dataset that received miraclib treatment, 
        and summarizing the baseline subgroup.
        """
    )

    st.subheader("Analysis workflow")

    wf1, wf2, wf3, wf4, wf5 = st.columns(5)
    wf1.markdown("### 1\n**Raw CSV**\n\nClinical sample metadata and five immune cell counts per sample")
    wf2.markdown("### 2\n**SQLite database**\n\nNormalized tables for projects, subjects, samples, populations, and counts")
    wf3.markdown("### 3\n**Frequency Summary**\n\nRelative frequency table for each population within each sample")
    wf4.markdown("### 4\n**Response Comparison**\n\nResponder vs non-responder comparison in melanoma PBMC samples treated with miraclib")
    wf5.markdown("### 5\n**Baseline Subset**\n\nBaseline melanoma PBMC subgroup summary and B-cell result")

    st.markdown("---")

    st.subheader("What each tab shows")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown(
            """
            **Overview**  
            Explains the workflow, generated outputs, and how to use the dashboard.

            **Part 2: Frequency Summary**  
            Shows the relative frequency of each immune cell population in each sample.
            """
        )
    with col2:
        st.markdown(
            """
            **Part 3: Statistical Analysis**  
            Compares responders and non-responders in melanoma PBMC samples treated with
            miraclib using a repeated-measures GEE model and interactive boxplots.

            **Part 4: Baseline Subset**  
            Summarizes the baseline melanoma PBMC miraclib subset, including project counts,
            response counts, sex counts, and the average B-cell count in melanoma male responders.
            """
        )

    st.markdown("---")

    st.subheader("Project snapshot")

    snap1, snap2, snap3, snap4 = st.columns(4)
    snap1.metric("Total samples", "10,500")
    snap2.metric("Immune populations", "5")
    snap3.metric("Part 3 samples", "1,968")
    snap4.metric("Part 4 baseline samples", "656")

    with st.expander("Database schema"):
        st.image("assets/schema.svg", caption="Normalized SQLite schema", use_container_width=True)

    st.markdown("---")

    st.subheader("Generated output status")

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
            key="part2_max_rows",
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

    if part3_stats_df is None or part3_filtered_df is None:
        st.warning("Part 3 output files were not found. Run the pipeline first.")
    else:
        st.markdown(
            """
            This analysis focuses on melanoma PBMC samples from subjects treated with miraclib.
            Response groups are compared using a binomial generalized estimating equation (GEE)
            model with repeated measures at timepoints 0, 7, and 14.
            """
        )

        responder_subjects = 0
        non_responder_subjects = 0
        unique_samples = 0

        if "subject_uid" in part3_filtered_df.columns:
            responder_subjects = (
                part3_filtered_df.loc[part3_filtered_df["response"] == "yes", "subject_uid"]
                .drop_duplicates()
                .shape[0]
            )
            non_responder_subjects = (
                part3_filtered_df.loc[part3_filtered_df["response"] == "no", "subject_uid"]
                .drop_duplicates()
                .shape[0]
            )

        if "sample" in part3_filtered_df.columns:
            unique_samples = part3_filtered_df["sample"].nunique()

        significant_count = 0
        min_fdr = None

        if "significant_fdr_0_05" in part3_stats_df.columns:
            significant_count = int(part3_stats_df["significant_fdr_0_05"].sum())

        if "joint_response_time_fdr_bh" in part3_stats_df.columns:
            min_fdr = part3_stats_df["joint_response_time_fdr_bh"].min()

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Populations tested", f"{len(part3_stats_df):,}")
        col2.metric("Unique samples", f"{unique_samples:,}")
        col3.metric("Responder subjects", f"{responder_subjects:,}")
        col4.metric("Non-responder subjects", f"{non_responder_subjects:,}")

        col5, col6 = st.columns(2)
        col5.metric("FDR-significant populations", f"{significant_count}")
        col6.metric(
            "Smallest FDR",
            "N/A" if min_fdr is None or pd.isna(min_fdr) else f"{min_fdr:.6f}",
        )

        st.subheader("Statistical results")
        st.dataframe(
            format_part3_stats_table(part3_stats_df),
            use_container_width=True,
            hide_index=True,
        )

        st.subheader("Interpretation")
        st.info(build_part3_interpretation(part3_stats_df))

        st.subheader("Responder vs non-responder boxplots")
        st.plotly_chart(
            build_part3_interactive_boxplot(part3_filtered_df),
            use_container_width=True,
        )

        if part3_summary_text is not None:
            with st.expander("Full GEE model summaries"):
                st.code(part3_summary_text, language="text")

with tab_part4:
    st.header("Part 4: Baseline Subset")

    if (
        part4_baseline_df is None
        or part4_project_df is None
        or part4_response_df is None
        or part4_sex_df is None
    ):
        st.warning("Part 4 output files not found. Run the pipeline first.")
    else:
        st.markdown(
            """
            This section focuses on melanoma PBMC samples collected at baseline
            (`time_from_treatment_start = 0`) from subjects treated with miraclib.
            """
        )

        baseline_samples = len(part4_baseline_df)
        unique_subjects = part4_baseline_df["subject"].nunique() if "subject" in part4_baseline_df.columns else 0
        n_projects = part4_baseline_df["project"].nunique() if "project" in part4_baseline_df.columns else 0

        male_responder_df = part4_baseline_df[
            (part4_baseline_df["sex"] == "M") & (part4_baseline_df["response"] == "yes")
        ].copy()

        avg_b_cell = None
        if not male_responder_df.empty and "b_cell_count" in male_responder_df.columns:
            avg_b_cell = male_responder_df["b_cell_count"].mean()

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Baseline samples", f"{baseline_samples:,}")
        col2.metric("Unique subjects", f"{unique_subjects:,}")
        col3.metric("Projects represented", f"{n_projects}")
        col4.metric(
            "Avg B-cell count (male responders)",
            "N/A" if avg_b_cell is None or pd.isna(avg_b_cell) else f"{avg_b_cell:.2f}",
        )

        st.subheader("Key Part 4 result")
        if avg_b_cell is not None and not pd.isna(avg_b_cell):
            st.success(
                "Considering melanoma males, the average number of B cells for responders at time = 0 is "
                f"{avg_b_cell:.2f}."
            )
        else:
            st.info("Average B-cell count for melanoma male responders could not be calculated.")

        table_col1, table_col2, table_col3 = st.columns(3)

        with table_col1:
            st.subheader("Samples per project")
            st.dataframe(part4_project_df, use_container_width=True, hide_index=True)

        with table_col2:
            st.subheader("Subjects by response")
            st.dataframe(part4_response_df, use_container_width=True, hide_index=True)

        with table_col3:
            st.subheader("Subjects by sex")
            st.dataframe(part4_sex_df, use_container_width=True, hide_index=True)

        st.subheader("Baseline subset table")

        project_options = sorted(part4_baseline_df["project"].dropna().unique().tolist())
        response_options = sorted(part4_baseline_df["response"].dropna().unique().tolist())
        sex_options = sorted(part4_baseline_df["sex"].dropna().unique().tolist())

        filter_col1, filter_col2, filter_col3 = st.columns(3)

        with filter_col1:
            selected_projects = st.multiselect(
                "Project",
                options=project_options,
                default=project_options,
            )

        with filter_col2:
            selected_responses = st.multiselect(
                "Response",
                options=response_options,
                default=response_options,
            )

        with filter_col3:
            selected_sexes = st.multiselect(
                "Sex",
                options=sex_options,
                default=sex_options,
            )

        max_rows_part4 = st.slider(
            "Maximum baseline rows to display",
            min_value=10,
            max_value=500,
            value=100,
            step=10,
            key="part4_max_rows",
        )

        filtered_part4_df = part4_baseline_df.copy()

        if selected_projects:
            filtered_part4_df = filtered_part4_df[
                filtered_part4_df["project"].isin(selected_projects)
            ]

        if selected_responses:
            filtered_part4_df = filtered_part4_df[
                filtered_part4_df["response"].isin(selected_responses)
            ]

        if selected_sexes:
            filtered_part4_df = filtered_part4_df[
                filtered_part4_df["sex"].isin(selected_sexes)
            ]

        st.write(
            f"Showing {min(len(filtered_part4_df), max_rows_part4):,} of {len(filtered_part4_df):,} rows"
        )

        st.dataframe(
            filtered_part4_df.head(max_rows_part4),
            use_container_width=True,
            hide_index=True,
        )

        if part4_summary_text is not None:
            with st.expander("Part 4 summary text"):
                st.code(part4_summary_text, language="text")
