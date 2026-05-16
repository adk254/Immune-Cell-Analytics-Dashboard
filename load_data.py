"""
load_data.py

Initializes a SQLite database and loads immune cell count data from cell-count.csv.

Run from the repository root with:

    python load_data.py
"""

from pathlib import Path
import sqlite3

import pandas as pd


ROOT_DIR = Path(__file__).resolve().parent
CSV_PATH = ROOT_DIR / "cell-count.csv"
DB_PATH = ROOT_DIR / "teiko_trial.db"

CELL_POPULATIONS = [
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
]


def create_schema(conn: sqlite3.Connection) -> None:
    """Create database schema for sample metadata and immune cell counts."""

    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS cell_counts;")
    cursor.execute("DROP TABLE IF EXISTS samples;")

    cursor.execute(
        """
        CREATE TABLE samples (
            sample TEXT PRIMARY KEY,
            project TEXT NOT NULL,
            subject TEXT NOT NULL,
            condition TEXT NOT NULL,
            age INTEGER,
            sex TEXT,
            treatment TEXT,
            response TEXT,
            sample_type TEXT,
            time_from_treatment_start INTEGER
        );
        """
    )

    cursor.execute(
        """
        CREATE TABLE cell_counts (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            sample TEXT NOT NULL,
            population TEXT NOT NULL,
            count INTEGER NOT NULL,
            FOREIGN KEY (sample) REFERENCES samples(sample)
        );
        """
    )

    cursor.execute(
        """
        CREATE INDEX idx_cell_counts_sample
        ON cell_counts(sample);
        """
    )

    cursor.execute(
        """
        CREATE INDEX idx_cell_counts_population
        ON cell_counts(population);
        """
    )

    cursor.execute(
        """
        CREATE INDEX idx_samples_trial_filters
        ON samples(condition, treatment, sample_type, response, time_from_treatment_start);
        """
    )

    conn.commit()


def load_data(conn: sqlite3.Connection, csv_path: Path) -> None:
    """Load cell-count.csv into the SQLite database."""

    if not csv_path.exists():
        raise FileNotFoundError(
            f"Could not find {csv_path}. Make sure cell-count.csv is in the repository root."
        )

    df = pd.read_csv(csv_path)

    required_columns = {
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
        *CELL_POPULATIONS,
    }

    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    sample_columns = [
        "sample",
        "project",
        "subject",
        "condition",
        "age",
        "sex",
        "treatment",
        "response",
        "sample_type",
        "time_from_treatment_start",
    ]

    samples_df = df[sample_columns].drop_duplicates(subset=["sample"])

    if len(samples_df) != len(df):
        print(
            "Warning: duplicate sample IDs were found. "
            "Only one metadata row per sample was loaded into the samples table."
        )

    samples_df.to_sql("samples", conn, if_exists="append", index=False)

    counts_df = df.melt(
        id_vars=["sample"],
        value_vars=CELL_POPULATIONS,
        var_name="population",
        value_name="count",
    )

    counts_df.to_sql("cell_counts", conn, if_exists="append", index=False)

    conn.commit()


def main() -> None:
    """Create database and load data."""

    conn = sqlite3.connect(DB_PATH)

    try:
        create_schema(conn)
        load_data(conn, CSV_PATH)

        sample_count = pd.read_sql_query("SELECT COUNT(*) AS n FROM samples;", conn)[
            "n"
        ].iloc[0]
        cell_count_rows = pd.read_sql_query(
            "SELECT COUNT(*) AS n FROM cell_counts;", conn
        )["n"].iloc[0]

        print(f"Database created: {DB_PATH}")
        print(f"Loaded {sample_count:,} samples.")
        print(f"Loaded {cell_count_rows:,} cell-count records.")

    finally:
        conn.close()


if __name__ == "__main__":
    main()
