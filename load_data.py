import sqlite3
from pathlib import Path

import pandas as pd

CSV_PATH = Path("cell-count.csv")
DB_PATH = Path("trial_data.db")

POPULATION_COLUMNS = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]

REQUIRED_COLUMNS = [
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    *POPULATION_COLUMNS,
]


def validate_dataframe(df: pd.DataFrame) -> None:
    missing_columns = [col for col in REQUIRED_COLUMNS if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    if df["sample"].duplicated().any():
        duplicate_samples = df.loc[df["sample"].duplicated(), "sample"].tolist()[:10]
        raise ValueError(f"Duplicate sample IDs found, example(s): {duplicate_samples}")

    metadata_columns = ["condition", "age", "sex", "treatment", "response"]
    inconsistent = []

    grouped = df.groupby(["project", "subject"], dropna=False)

    for col in metadata_columns:
        n_unique = grouped[col].nunique(dropna=False)
        bad_groups = n_unique[n_unique > 1]
        if not bad_groups.empty:
            inconsistent.append(col)

    if inconsistent:
        raise ValueError(
            "Subject-level metadata is inconsistent within (project, subject) for: "
            + ", ".join(inconsistent)
        )

    for col in POPULATION_COLUMNS:
        if not pd.api.types.is_numeric_dtype(df[col]):
            raise ValueError(f"Column '{col}' must be numeric.")
        if (df[col] < 0).any():
            raise ValueError(f"Column '{col}' contains negative values.")

    if not pd.api.types.is_numeric_dtype(df["time_from_treatment_start"]):
        raise ValueError("Column 'time_from_treatment_start' must be numeric.")

    # Make sure time values are whole numbers
    time_values = pd.to_numeric(df["time_from_treatment_start"], errors="coerce")
    if time_values.isna().any():
        raise ValueError("Column 'time_from_treatment_start' contains non-numeric values.")

    if not (time_values == time_values.astype(int)).all():
        raise ValueError("Column 'time_from_treatment_start' must contain integer values.")


def initialize_database(conn: sqlite3.Connection) -> None:
    conn.execute("PRAGMA foreign_keys = ON;")

    conn.executescript(
        """
        DROP TABLE IF EXISTS cell_counts;
        DROP TABLE IF EXISTS populations;
        DROP TABLE IF EXISTS samples;
        DROP TABLE IF EXISTS subjects;
        DROP TABLE IF EXISTS projects;

        CREATE TABLE projects (
            project_id INTEGER PRIMARY KEY,
            project_code TEXT NOT NULL UNIQUE
        );

        CREATE TABLE subjects (
            subject_id INTEGER PRIMARY KEY,
            project_id INTEGER NOT NULL,
            subject_code TEXT NOT NULL,
            condition TEXT NOT NULL,
            age INTEGER NOT NULL,
            sex TEXT NOT NULL,
            treatment TEXT NOT NULL,
            response TEXT NULL,
            UNIQUE(project_id, subject_code),
            FOREIGN KEY (project_id) REFERENCES projects(project_id)
        );

        CREATE TABLE samples (
            sample_id INTEGER PRIMARY KEY,
            sample_code TEXT NOT NULL UNIQUE,
            subject_id INTEGER NOT NULL,
            sample_type TEXT NOT NULL,
            time_from_treatment_start INTEGER NOT NULL,
            FOREIGN KEY (subject_id) REFERENCES subjects(subject_id)
        );

        CREATE TABLE populations (
            population_id INTEGER PRIMARY KEY,
            population_name TEXT NOT NULL UNIQUE
        );

        CREATE TABLE cell_counts (
            sample_id INTEGER NOT NULL,
            population_id INTEGER NOT NULL,
            count INTEGER NOT NULL,
            PRIMARY KEY (sample_id, population_id),
            FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
            FOREIGN KEY (population_id) REFERENCES populations(population_id)
        );
        """
    )


def load_projects(conn: sqlite3.Connection, df: pd.DataFrame) -> dict[str, int]:
    project_codes = sorted(df["project"].drop_duplicates().tolist())
    conn.executemany(
        "INSERT INTO projects (project_code) VALUES (?);",
        [(code,) for code in project_codes],
    )

    rows = conn.execute("SELECT project_id, project_code FROM projects;").fetchall()
    return {project_code: project_id for project_id, project_code in rows}


def load_subjects(
    conn: sqlite3.Connection, df: pd.DataFrame, project_map: dict[str, int]
) -> dict[tuple[str, str], int]:
    subject_df = (
        df[
            [
                "project",
                "subject",
                "condition",
                "age",
                "sex",
                "treatment",
                "response",
            ]
        ]
        .drop_duplicates()
        .copy()
    )

    subject_records = []
    for row in subject_df.itertuples(index=False):
        subject_records.append(
            (
                project_map[row.project],
                row.subject,
                row.condition,
                int(row.age),
                row.sex,
                row.treatment,
                None if pd.isna(row.response) else row.response,
            )
        )

    conn.executemany(
        """
        INSERT INTO subjects (
            project_id, subject_code, condition, age, sex, treatment, response
        ) VALUES (?, ?, ?, ?, ?, ?, ?);
        """,
        subject_records,
    )

    rows = conn.execute(
        """
        SELECT s.subject_id, p.project_code, s.subject_code
        FROM subjects s
        JOIN projects p ON s.project_id = p.project_id;
        """
    ).fetchall()

    return {(project_code, subject_code): subject_id for subject_id, project_code, subject_code in rows}


def load_populations(conn: sqlite3.Connection) -> dict[str, int]:
    conn.executemany(
        "INSERT INTO populations (population_name) VALUES (?);",
        [(name,) for name in POPULATION_COLUMNS],
    )

    rows = conn.execute("SELECT population_id, population_name FROM populations;").fetchall()
    return {population_name: population_id for population_id, population_name in rows}


def load_samples(
    conn: sqlite3.Connection,
    df: pd.DataFrame,
    subject_map: dict[tuple[str, str], int],
) -> dict[str, int]:
    sample_records = []
    for row in df[
        ["project", "subject", "sample", "sample_type", "time_from_treatment_start"]
    ].itertuples(index=False):
        subject_id = subject_map[(row.project, row.subject)]
        sample_records.append(
            (
                row.sample,
                subject_id,
                row.sample_type,
                int(row.time_from_treatment_start),
            )
        )

    conn.executemany(
        """
        INSERT INTO samples (
            sample_code, subject_id, sample_type, time_from_treatment_start
        ) VALUES (?, ?, ?, ?);
        """,
        sample_records,
    )

    rows = conn.execute("SELECT sample_id, sample_code FROM samples;").fetchall()
    return {sample_code: sample_id for sample_id, sample_code in rows}


def load_cell_counts(
    conn: sqlite3.Connection,
    df: pd.DataFrame,
    sample_map: dict[str, int],
    population_map: dict[str, int],
) -> None:
    count_records = []

    for row in df[["sample", *POPULATION_COLUMNS]].itertuples(index=False):
        sample_code = row[0]
        sample_id = sample_map[sample_code]

        for population_name, count_value in zip(POPULATION_COLUMNS, row[1:]):
            count_records.append(
                (
                    sample_id,
                    population_map[population_name],
                    int(count_value),
                )
            )

    conn.executemany(
        """
        INSERT INTO cell_counts (sample_id, population_id, count)
        VALUES (?, ?, ?);
        """,
        count_records,
    )


def main() -> None:
    if not CSV_PATH.exists():
        raise FileNotFoundError(f"Could not find input file: {CSV_PATH}")

    df = pd.read_csv(CSV_PATH)
    validate_dataframe(df)

    with sqlite3.connect(DB_PATH) as conn:
        initialize_database(conn)

        project_map = load_projects(conn, df)
        subject_map = load_subjects(conn, df, project_map)
        population_map = load_populations(conn)
        sample_map = load_samples(conn, df, subject_map)
        load_cell_counts(conn, df, sample_map, population_map)

        conn.commit()

        n_projects = conn.execute("SELECT COUNT(*) FROM projects;").fetchone()[0]
        n_subjects = conn.execute("SELECT COUNT(*) FROM subjects;").fetchone()[0]
        n_samples = conn.execute("SELECT COUNT(*) FROM samples;").fetchone()[0]
        n_populations = conn.execute("SELECT COUNT(*) FROM populations;").fetchone()[0]
        n_counts = conn.execute("SELECT COUNT(*) FROM cell_counts;").fetchone()[0]

    print(f"Database created: {DB_PATH}")
    print(f"Projects loaded: {n_projects}")
    print(f"Subjects loaded: {n_subjects}")
    print(f"Samples loaded: {n_samples}")
    print(f"Populations loaded: {n_populations}")
    print(f"Cell count rows loaded: {n_counts}")


if __name__ == "__main__":
    main()
