import sqlite3
from pathlib import Path

import pandas as pd

DB_PATH = Path("trial_data.db")
OUTPUT_DIR = Path("outputs")
OUTPUT_PATH = OUTPUT_DIR / "sample_population_frequencies.csv"


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


def build_frequency_summary(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    df["total_count"] = df.groupby("sample")["count"].transform("sum")
    df["percentage"] = (df["count"] / df["total_count"]) * 100

    summary_df = df[["sample", "total_count", "population", "count", "percentage"]].copy()

    summary_df = summary_df.sort_values(["sample", "population"]).reset_index(drop=True)
    return summary_df


def main() -> None:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Could not find database file: {DB_PATH}. Run python load_data.py first."
        )

    OUTPUT_DIR.mkdir(exist_ok=True)

    with sqlite3.connect(DB_PATH) as conn:
        counts_df = fetch_cell_count_data(conn)

    summary_df = build_frequency_summary(counts_df)
    summary_df.to_csv(OUTPUT_PATH, index=False)

    print(f"Summary table written to: {OUTPUT_PATH}")
    print(f"Rows written: {len(summary_df)}")
    print("\nFirst 10 rows:")
    print(summary_df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
