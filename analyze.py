"""
analyze.py

Runs the full analysis pipeline after the SQLite database has been created.

Run from the repository root with:

    python analyze.py
"""

from src.analysis import (
    get_baseline_melanoma_miraclib_pbmc_samples,
    get_baseline_project_counts,
    get_baseline_response_subject_counts,
    get_baseline_sex_subject_counts,
    get_cell_frequency_table,
    get_male_responder_baseline_bcell_average,
    get_melanoma_miraclib_pbmc_frequencies,
)
from src.db import OUTPUTS_DIR, ensure_outputs_dir
from src.plots import create_response_boxplot
from src.stats import compare_response_groups


def main() -> None:
    """Generate all required output tables, statistics, and plots."""

    ensure_outputs_dir()

    print("Generating cell population frequency table...")
    frequency_df = get_cell_frequency_table()
    frequency_output = OUTPUTS_DIR / "cell_frequencies.csv"
    frequency_df.to_csv(frequency_output, index=False)
    print(f"Saved {frequency_output}")

    print("Filtering melanoma PBMC samples treated with miraclib...")
    response_frequency_df = get_melanoma_miraclib_pbmc_frequencies()
    response_frequency_output = OUTPUTS_DIR / "melanoma_miraclib_pbmc_frequencies.csv"
    response_frequency_df.to_csv(response_frequency_output, index=False)
    print(f"Saved {response_frequency_output}")

    print("Running responder vs non-responder statistical tests...")
    stats_df = compare_response_groups(response_frequency_df)
    stats_output = OUTPUTS_DIR / "responder_vs_nonresponder_stats.csv"
    stats_df.to_csv(stats_output, index=False)
    print(f"Saved {stats_output}")

    print("Creating responder vs non-responder boxplot...")
    boxplot_output = OUTPUTS_DIR / "responder_vs_nonresponder_boxplot.html"
    create_response_boxplot(response_frequency_df, boxplot_output)
    print(f"Saved {boxplot_output}")

    print("Generating baseline subset outputs...")
    baseline_samples_df = get_baseline_melanoma_miraclib_pbmc_samples()
    baseline_samples_output = OUTPUTS_DIR / "baseline_melanoma_miraclib_pbmc_samples.csv"
    baseline_samples_df.to_csv(baseline_samples_output, index=False)
    print(f"Saved {baseline_samples_output}")

    project_counts_df = get_baseline_project_counts()
    project_counts_output = OUTPUTS_DIR / "baseline_project_counts.csv"
    project_counts_df.to_csv(project_counts_output, index=False)
    print(f"Saved {project_counts_output}")

    response_counts_df = get_baseline_response_subject_counts()
    response_counts_output = OUTPUTS_DIR / "baseline_response_subject_counts.csv"
    response_counts_df.to_csv(response_counts_output, index=False)
    print(f"Saved {response_counts_output}")

    sex_counts_df = get_baseline_sex_subject_counts()
    sex_counts_output = OUTPUTS_DIR / "baseline_sex_subject_counts.csv"
    sex_counts_df.to_csv(sex_counts_output, index=False)
    print(f"Saved {sex_counts_output}")

    avg_bcell = get_male_responder_baseline_bcell_average()
    avg_bcell_output = OUTPUTS_DIR / "male_responder_baseline_bcell_average.txt"

    with open(avg_bcell_output, "w", encoding="utf-8") as file:
        file.write(f"{avg_bcell:.2f}\n")

    print(f"Saved {avg_bcell_output}")
    print(f"Average B-cell count for male melanoma responders at baseline: {avg_bcell:.2f}")

    print("Analysis pipeline completed successfully.")


if __name__ == "__main__":
    main()