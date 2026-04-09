import sqlite3
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind
from statsmodels.stats.multitest import multipletests

DB_PATH = Path("trial_data.db")
OUTPUT_DIR = Path("outputs")

PART2_OUTPUT = OUTPUT_DIR / "sample_population_frequencies.csv"
PART3_FILTERED_OUTPUT = OUTPUT_DIR / "part3_filtered_frequencies.csv"
PART3_WELCH_OUTPUT = OUTPUT_DIR / "part3_population_statistics_welch.csv"
PART3_MWU_OUTPUT = OUTPUT_DIR / "part3_population_statistics_mannwhitney.csv"
PART3_PLOT_OUTPUT = OUTPUT_DIR / "part3_boxplots.png"

PART4_BASELINE_OUTPUT = OUTPUT_DIR / "part4_baseline_samples.csv"
PART4_PROJECT_OUTPUT = OUTPUT_DIR / "part4_samples_per_project.csv"
PART4_RESPONSE_OUTPUT = OUTPUT_DIR / "part4_subjects_by_response.csv"
PART4_SEX_OUTPUT = OUTPUT_DIR / "part4_subjects_by_sex.csv"
PART4_SUMMARY_OUTPUT = OUTPUT_DIR / "part4_summary.txt"


def fetch_cell_count_data(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
    SELECT
        s.sample_code AS sample,
        p.population_name AS population,
        cc.count
    FROM cell_counts cc
    JOIN samples s
        ON cc.sample_id = s.sample_id
    JOIN populations p
        ON cc.population_id = p.population_id
    ORDER BY s.sample_code, p.population_name;
    """
    return pd.read_sql_query(query, conn)


def fetch_sample_metadata(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
    SELECT
        s.sample_code AS sample,
        pr.project_code AS project,
        sub.subject_code AS subject,
        sub.condition,
        sub.age,
        sub.sex,
        sub.treatment,
        sub.response,
        s.sample_type,
        s.time_from_treatment_start
    FROM samples s
    JOIN subjects sub
        ON s.subject_id = sub.subject_id
    JOIN projects pr
        ON sub.project_id = pr.project_id
    ORDER BY s.sample_code;
    """
    return pd.read_sql_query(query, conn)


def fetch_b_cell_counts(conn: sqlite3.Connection) -> pd.DataFrame:
    query = """
    SELECT
        s.sample_code AS sample,
        cc.count AS b_cell_count
    FROM cell_counts cc
    JOIN samples s
        ON cc.sample_id = s.sample_id
    JOIN populations p
        ON cc.population_id = p.population_id
    WHERE p.population_name = 'b_cell'
    ORDER BY s.sample_code;
    """
    return pd.read_sql_query(query, conn)


def build_frequency_summary(counts_df: pd.DataFrame) -> pd.DataFrame:
    summary_df = counts_df.copy()
    summary_df["total_count"] = summary_df.groupby("sample")["count"].transform("sum")
    summary_df["percentage"] = (summary_df["count"] / summary_df["total_count"]) * 100

    summary_df = summary_df[
        ["sample", "total_count", "population", "count", "percentage"]
    ].copy()

    summary_df = summary_df.sort_values(["sample", "population"]).reset_index(drop=True)
    return summary_df


def filter_part3_data(
    summary_df: pd.DataFrame, metadata_df: pd.DataFrame
) -> pd.DataFrame:
    merged_df = summary_df.merge(metadata_df, on="sample", how="left")

    filtered_df = merged_df[
        (merged_df["condition"] == "melanoma")
        & (merged_df["treatment"] == "miraclib")
        & (merged_df["sample_type"] == "PBMC")
        & (merged_df["response"].isin(["yes", "no"]))
    ].copy()

    filtered_df = filtered_df.sort_values(
        ["population", "response", "sample"]
    ).reset_index(drop=True)
    return filtered_df


def summarize_groups(
    responder_vals: pd.Series, non_responder_vals: pd.Series
) -> dict:
    return {
        "responder_n": len(responder_vals),
        "non_responder_n": len(non_responder_vals),
        "responder_mean_percentage": responder_vals.mean(),
        "non_responder_mean_percentage": non_responder_vals.mean(),
        "responder_median_percentage": responder_vals.median(),
        "non_responder_median_percentage": non_responder_vals.median(),
        "mean_difference_percentage": responder_vals.mean() - non_responder_vals.mean(),
        "median_difference_percentage": responder_vals.median() - non_responder_vals.median(),
    }


def run_part3_welch_statistics(filtered_df: pd.DataFrame) -> pd.DataFrame:
    results = []
    populations = sorted(filtered_df["population"].unique().tolist())

    for population in populations:
        pop_df = filtered_df[filtered_df["population"] == population]

        responder_vals = pop_df.loc[pop_df["response"] == "yes", "percentage"]
        non_responder_vals = pop_df.loc[pop_df["response"] == "no", "percentage"]

        statistic, p_value = ttest_ind(
            responder_vals,
            non_responder_vals,
            equal_var=False,
            nan_policy="raise",
        )

        results.append(
            {
                "population": population,
                **summarize_groups(responder_vals, non_responder_vals),
                "test": "welch_t_test",
                "test_statistic": statistic,
                "p_value": p_value,
            }
        )

    stats_df = pd.DataFrame(results)
    rejected, fdr_values, _, _ = multipletests(stats_df["p_value"], method="fdr_bh")
    stats_df["fdr_bh"] = fdr_values
    stats_df["significant_fdr_0_05"] = rejected

    stats_df = stats_df.sort_values("p_value").reset_index(drop=True)
    return stats_df


def run_part3_mannwhitney_statistics(filtered_df: pd.DataFrame) -> pd.DataFrame:
    results = []
    populations = sorted(filtered_df["population"].unique().tolist())

    for population in populations:
        pop_df = filtered_df[filtered_df["population"] == population]

        responder_vals = pop_df.loc[pop_df["response"] == "yes", "percentage"]
        non_responder_vals = pop_df.loc[pop_df["response"] == "no", "percentage"]

        statistic, p_value = mannwhitneyu(
            responder_vals,
            non_responder_vals,
            alternative="two-sided",
        )

        results.append(
            {
                "population": population,
                **summarize_groups(responder_vals, non_responder_vals),
                "test": "mann_whitney_u",
                "test_statistic": statistic,
                "p_value": p_value,
            }
        )

    stats_df = pd.DataFrame(results)
    rejected, fdr_values, _, _ = multipletests(stats_df["p_value"], method="fdr_bh")
    stats_df["fdr_bh"] = fdr_values
    stats_df["significant_fdr_0_05"] = rejected

    stats_df = stats_df.sort_values("p_value").reset_index(drop=True)
    return stats_df


def make_part3_boxplots(filtered_df: pd.DataFrame, output_path: Path) -> None:
    populations = sorted(filtered_df["population"].unique().tolist())
    fig, axes = plt.subplots(1, len(populations), figsize=(20, 5), sharey=True)

    if len(populations) == 1:
        axes = [axes]

    for ax, population in zip(axes, populations):
        pop_df = filtered_df[filtered_df["population"] == population]

        non_responder_vals = pop_df.loc[pop_df["response"] == "no", "percentage"]
        responder_vals = pop_df.loc[pop_df["response"] == "yes", "percentage"]

        ax.boxplot([non_responder_vals, responder_vals], tick_labels=["No", "Yes"])
        ax.set_title(population)
        ax.set_xlabel("Response")

    axes[0].set_ylabel("Relative frequency (%)")
    fig.suptitle(
        "Melanoma PBMC samples treated with miraclib: responder vs non-responder"
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_part4_baseline_subset(
    metadata_df: pd.DataFrame, b_cell_df: pd.DataFrame
) -> pd.DataFrame:
    baseline_df = metadata_df.merge(b_cell_df, on="sample", how="left")

    baseline_df = baseline_df[
        (baseline_df["condition"] == "melanoma")
        & (baseline_df["treatment"] == "miraclib")
        & (baseline_df["sample_type"] == "PBMC")
        & (baseline_df["time_from_treatment_start"] == 0)
    ].copy()

    baseline_df = baseline_df.sort_values(["project", "subject", "sample"]).reset_index(drop=True)
    return baseline_df


def run_part4_outputs(baseline_df: pd.DataFrame) -> dict:
    samples_per_project = (
        baseline_df.groupby("project", as_index=False)
        .agg(n_samples=("sample", "count"))
        .sort_values("project")
        .reset_index(drop=True)
    )

    subjects_by_response = (
        baseline_df[["project", "subject", "response"]]
        .drop_duplicates()
        .groupby("response", as_index=False)
        .agg(n_subjects=("subject", "count"))
        .sort_values("response")
        .reset_index(drop=True)
    )

    subjects_by_sex = (
        baseline_df[["project", "subject", "sex"]]
        .drop_duplicates()
        .groupby("sex", as_index=False)
        .agg(n_subjects=("subject", "count"))
        .sort_values("sex")
        .reset_index(drop=True)
    )

    male_responder_df = baseline_df[
        (baseline_df["sex"] == "M") & (baseline_df["response"] == "yes")
    ].copy()

    average_b_cell_count = round(male_responder_df["b_cell_count"].mean(), 2)

    return {
        "samples_per_project": samples_per_project,
        "subjects_by_response": subjects_by_response,
        "subjects_by_sex": subjects_by_sex,
        "average_b_cell_count_male_responders": average_b_cell_count,
    }


def main() -> None:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Could not find database file: {DB_PATH}. Run python load_data.py first."
        )

    OUTPUT_DIR.mkdir(exist_ok=True)

    with sqlite3.connect(DB_PATH) as conn:
        counts_df = fetch_cell_count_data(conn)
        metadata_df = fetch_sample_metadata(conn)
        b_cell_df = fetch_b_cell_counts(conn)

    # Part 2
    summary_df = build_frequency_summary(counts_df)
    summary_df.to_csv(PART2_OUTPUT, index=False)

    # Part 3
    filtered_df = filter_part3_data(summary_df, metadata_df)
    filtered_df.to_csv(PART3_FILTERED_OUTPUT, index=False)

    welch_df = run_part3_welch_statistics(filtered_df)
    welch_df.to_csv(PART3_WELCH_OUTPUT, index=False)

    mwu_df = run_part3_mannwhitney_statistics(filtered_df)
    mwu_df.to_csv(PART3_MWU_OUTPUT, index=False)

    make_part3_boxplots(filtered_df, PART3_PLOT_OUTPUT)

    # Part 4
    baseline_df = build_part4_baseline_subset(metadata_df, b_cell_df)
    baseline_df.to_csv(PART4_BASELINE_OUTPUT, index=False)

    part4_results = run_part4_outputs(baseline_df)

    part4_results["samples_per_project"].to_csv(PART4_PROJECT_OUTPUT, index=False)
    part4_results["subjects_by_response"].to_csv(PART4_RESPONSE_OUTPUT, index=False)
    part4_results["subjects_by_sex"].to_csv(PART4_SEX_OUTPUT, index=False)

    with open(PART4_SUMMARY_OUTPUT, "w", encoding="utf-8") as f:
        f.write("Part 4 baseline subset summary\n")
        f.write("==============================\n")
        f.write(f"Baseline melanoma PBMC miraclib samples: {len(baseline_df)}\n")
        f.write(
            "Average B-cell count for melanoma male responders at time=0: "
            f"{part4_results['average_b_cell_count_male_responders']:.2f}\n"
        )

    unique_samples = (
        filtered_df[["sample", "response"]]
        .drop_duplicates()
        .groupby("response")
        .size()
        .to_dict()
    )

    print(f"Part 2 summary table written to: {PART2_OUTPUT}")
    print(f"Part 2 rows written: {len(summary_df)}")

    print(f"\nPart 3 filtered table written to: {PART3_FILTERED_OUTPUT}")
    print(f"Part 3 filtered rows written: {len(filtered_df)}")
    print(f"Unique filtered samples by response: {unique_samples}")

    print(f"\nPrimary analysis written to: {PART3_WELCH_OUTPUT}")
    print(welch_df.to_string(index=False))

    print(f"\nSensitivity analysis written to: {PART3_MWU_OUTPUT}")
    print(mwu_df.to_string(index=False))

    print(f"\nPart 3 boxplots written to: {PART3_PLOT_OUTPUT}")

    print(f"\nPart 4 baseline samples written to: {PART4_BASELINE_OUTPUT}")
    print(f"Baseline subset rows: {len(baseline_df)}")

    print(f"\nPart 4 samples per project written to: {PART4_PROJECT_OUTPUT}")
    print(part4_results["samples_per_project"].to_string(index=False))

    print(f"\nPart 4 subjects by response written to: {PART4_RESPONSE_OUTPUT}")
    print(part4_results["subjects_by_response"].to_string(index=False))

    print(f"\nPart 4 subjects by sex written to: {PART4_SEX_OUTPUT}")
    print(part4_results["subjects_by_sex"].to_string(index=False))

    print(f"\nPart 4 summary written to: {PART4_SUMMARY_OUTPUT}")
    print(
        "Average B-cell count for melanoma male responders at time=0: "
        f"{part4_results['average_b_cell_count_male_responders']:.2f}"
    )


if __name__ == "__main__":
    main()
