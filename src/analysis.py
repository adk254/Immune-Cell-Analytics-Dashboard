"""
Analysis functions for immune cell population summaries and clinical subsets.
"""

import pandas as pd

from src.db import read_query


def get_cell_frequency_table() -> pd.DataFrame:
    """
    Return relative frequency of each cell population within each sample.

    Output columns:
    - sample
    - total_count
    - population
    - count
    - percentage
    """

    query = """
        WITH sample_totals AS (
            SELECT
                sample,
                SUM(count) AS total_count
            FROM cell_counts
            GROUP BY sample
        )
        SELECT
            c.sample AS sample,
            t.total_count AS total_count,
            c.population AS population,
            c.count AS count,
            ROUND((CAST(c.count AS REAL) / t.total_count) * 100, 4) AS percentage
        FROM cell_counts c
        JOIN sample_totals t
            ON c.sample = t.sample
        ORDER BY c.sample, c.population;
    """

    return read_query(query)


def get_melanoma_miraclib_pbmc_frequencies() -> pd.DataFrame:
    """
    Return cell population frequencies for melanoma PBMC samples treated with miraclib.

    This table is used for responder vs non-responder statistical comparisons.
    """

    query = """
        WITH sample_totals AS (
            SELECT
                sample,
                SUM(count) AS total_count
            FROM cell_counts
            GROUP BY sample
        )
        SELECT
            s.project,
            s.subject,
            s.condition,
            s.age,
            s.sex,
            s.treatment,
            s.response,
            s.sample,
            s.sample_type,
            s.time_from_treatment_start,
            c.population,
            c.count,
            t.total_count,
            (CAST(c.count AS REAL) / t.total_count) * 100 AS percentage
        FROM samples s
        JOIN cell_counts c
            ON s.sample = c.sample
        JOIN sample_totals t
            ON s.sample = t.sample
        WHERE LOWER(s.condition) = 'melanoma'
            AND LOWER(s.treatment) = 'miraclib'
            AND UPPER(s.sample_type) = 'PBMC'
            AND LOWER(s.response) IN ('yes', 'no')
        ORDER BY s.sample, c.population;
    """

    return read_query(query)


def get_baseline_melanoma_miraclib_pbmc_samples() -> pd.DataFrame:
    """
    Return melanoma PBMC baseline samples from patients treated with miraclib.

    Baseline is defined as time_from_treatment_start = 0.
    """

    query = """
        SELECT
            project,
            subject,
            condition,
            age,
            sex,
            treatment,
            response,
            sample,
            sample_type,
            time_from_treatment_start
        FROM samples
        WHERE LOWER(condition) = 'melanoma'
            AND LOWER(treatment) = 'miraclib'
            AND UPPER(sample_type) = 'PBMC'
            AND time_from_treatment_start = 0
        ORDER BY project, subject, sample;
    """

    return read_query(query)


def get_baseline_project_counts() -> pd.DataFrame:
    """Return the number of baseline melanoma/miraclib/PBMC samples per project."""

    baseline_df = get_baseline_melanoma_miraclib_pbmc_samples()

    return (
        baseline_df.groupby("project", as_index=False)
        .agg(n_samples=("sample", "nunique"))
        .sort_values("project")
    )


def get_baseline_response_subject_counts() -> pd.DataFrame:
    """Return the number of unique baseline subjects by response status."""

    baseline_df = get_baseline_melanoma_miraclib_pbmc_samples()

    return (
        baseline_df.groupby("response", as_index=False)
        .agg(n_subjects=("subject", "nunique"))
        .sort_values("response")
    )


def get_baseline_sex_subject_counts() -> pd.DataFrame:
    """Return the number of unique baseline subjects by sex."""

    baseline_df = get_baseline_melanoma_miraclib_pbmc_samples()

    return (
        baseline_df.groupby("sex", as_index=False)
        .agg(n_subjects=("subject", "nunique"))
        .sort_values("sex")
    )


def get_male_responder_baseline_bcell_average() -> float:
    """
    Return average B-cell count for melanoma male responders at baseline
    who received miraclib and had PBMC samples.
    """

    query = """
        SELECT
            AVG(c.count) AS average_b_cell_count
        FROM samples s
        JOIN cell_counts c
            ON s.sample = c.sample
        WHERE LOWER(s.condition) = 'melanoma'
            AND LOWER(s.treatment) = 'miraclib'
            AND UPPER(s.sample_type) = 'PBMC'
            AND s.time_from_treatment_start = 0
            AND UPPER(s.sex) = 'M'
            AND LOWER(s.response) = 'yes'
            AND c.population = 'b_cell';
    """

    result = read_query(query)
    value = result["average_b_cell_count"].iloc[0]

    if pd.isna(value):
        return float("nan")

    return float(value)