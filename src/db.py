"""
Database utilities for the Teiko technical assessment.
"""

from pathlib import Path
import sqlite3

import pandas as pd


ROOT_DIR = Path(__file__).resolve().parents[1]
OUTPUTS_DIR = ROOT_DIR / "outputs"
DB_PATH = ROOT_DIR / "teiko_trial.db"


def get_connection() -> sqlite3.Connection:
    """Return a SQLite connection to the project database."""
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Database not found at {DB_PATH}. Run `python load_data.py` first."
        )

    return sqlite3.connect(DB_PATH)


def read_query(query: str, params: tuple | None = None) -> pd.DataFrame:
    """Run a SQL query and return the result as a pandas DataFrame."""
    with get_connection() as conn:
        return pd.read_sql_query(query, conn, params=params)


def ensure_outputs_dir() -> None:
    """Create the outputs directory if it does not already exist."""
    OUTPUTS_DIR.mkdir(exist_ok=True)
